/* Edoardo Caldarelli, Institut de Robòtica i Informàtica Industrial (CSIC-UPC)

code based on work by Adrià Luque Acera */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "general_functions.h"
#include "model_functions.h"
#include "mpc_functions.h"
#include <ros/ros.h>
#include "cartesian_msgs/CartesianCommand.h"
#include "geometry_msgs/Pose.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/Quaternion.h"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <mpc_pkg/TwoPoses.h>
#include <mpc_pkg/TwoPoints.h>
#include <mpc_pkg/OptiData.h>
#include <mpc_pkg/OptiDataQDMC.h>
#include <mpc_pkg/HorizonControls.h>
#include <mpc_vision/SOMstate.h>
#include <std_srvs/Empty.h>


using namespace std;


// GLOBAL VARIABLES
string datapath = "/home/robot/Desktop/ALuque/mpc_ros/src/mpc_node/data/";
int NTraj = 0;
int nCOM = 0;
int nSOM = 0;
int Hp = 0;
int Nss = 0;
int Hc = 0;
double Ts = 0.00;
Eigen::Vector3d tcp_offset_local(0.0, 0.0, 0.09);

// Needed to update Vision SOM
int MaxVSteps = 12;
double MaxVDiff = 0.03;
double Wv = 0.0;
Eigen::MatrixXd A_SOM;
Eigen::MatrixXd B_SOM;
Eigen::VectorXd f_SOM;
Eigen::VectorXd state_SOM;
Eigen::VectorXd SOM_nodes_ctrl(2);
Eigen::VectorXd SOM_coord_ctrl(6);
Eigen::MatrixXd uSOMhist(6,MaxVSteps);
int uhistID = 0;

// Needed to connect with optimizer
Eigen::MatrixXd uHp_SOM(Hp,6);
int usomhpID = 0;

// ----------------


// Linear model structure
typedef struct linmdl {
  int row;
  int col;
  double mass;
  double grav;
  double dt;
  Eigen::Vector3d stiffness;
  Eigen::Vector3d damping;
  double z_sum;
  Eigen::MatrixXd nodeInitial;
  Eigen::MatrixXd mat_x;
  Eigen::MatrixXd mat_y;
  Eigen::MatrixXd mat_z;
  Eigen::VectorXd coord_ctrl;
  Eigen::VectorXd coord_lc;
} LinMdl;



// ROS SUBSCRIBERS - CALLBACK FUNCTIONS
void somstateReceived(const mpc_vision::SOMstate& msg) {
	
	// Get timestamp and size
	double Vision_t = msg.header.stamp.toSec();
	int nMesh = msg.size;
	int Meshlength = nMesh*nMesh;
	
	// Obtain mesh from vision feedback
	Eigen::VectorXf Vision_SOMstatef = Eigen::VectorXf::Map(msg.states.data(), msg.states.size());
  Eigen::VectorXd Vision_SOMstate = Vision_SOMstatef.cast<double>();

	// Reduce to SOM mesh size
	Eigen::VectorXd VSOMpos(3*nSOM*nSOM);
  Eigen::VectorXd VSOMvel(3*nSOM*nSOM);
  Eigen::VectorXd VSOMstate(6*nSOM*nSOM);
  tie(VSOMpos, VSOMvel) = take_reduced_mesh(Vision_SOMstate.segment(0,3*Meshlength),
  																					Vision_SOMstate.segment(3*Meshlength,3*Meshlength),
  																					nMesh, nSOM);
  VSOMstate << VSOMpos, VSOMvel;
  																							
  //double Vision_dt = ros::Time::now().toSec() - Vision_t;
  ros::Duration Vision_dt = ros::Time::now() - msg.header.stamp;
  
  // Update the captured SOMstate to current time (dt/Ts steps)
  int update_steps = Vision_dt.toSec()/Ts;
  update_steps = max(update_steps, 0);
  if (update_steps > MaxVSteps) {
  	update_steps = MaxVSteps;
  	ROS_WARN_STREAM("Ignoring Vision data: Delay longer than maximum.");
  	return;
  }
  
  // Needed variables
  Eigen::VectorXd uSOM_Vti(6);
  Eigen::VectorXd ulin_Vti(6);
  Eigen::VectorXd urot_Vti(6);
  Eigen::MatrixXd ulin2_Vti(3,2);
  Eigen::MatrixXd urot2_Vti(3,2);
  Eigen::Vector3d Vcloth_x;
  Eigen::Vector3d Vcloth_y;
  Eigen::Vector3d Vcloth_z;
  Eigen::Matrix3d VRcloth;
  Eigen::MatrixXd pos_ini_VSOM(nSOM*nSOM,3);
  Eigen::MatrixXd vel_ini_VSOM(nSOM*nSOM,3);
  Eigen::MatrixXd pos_ini_VSOM_rot(nSOM*nSOM,3);
  Eigen::MatrixXd vel_ini_VSOM_rot(nSOM*nSOM,3);
  Eigen::MatrixXd pos_nxt_VSOM_rot(nSOM*nSOM,3);
  Eigen::MatrixXd vel_nxt_VSOM_rot(nSOM*nSOM,3);
  Eigen::MatrixXd pos_nxt_VSOM(nSOM*nSOM,3);
  Eigen::MatrixXd vel_nxt_VSOM(nSOM*nSOM,3);
  Eigen::VectorXd VSOMstate_rot(6*nSOM*nSOM);
	
	// Iterate until current time
	for (int Vti=-update_steps; Vti<0; Vti++) {
		
		// Real history ID
		int uhistIDi = uhistID + Vti;
		if (uhistIDi < 0) uhistIDi+=MaxVSteps;
		
		// Get uSOM from that instant, skip 0s
		uSOM_Vti = uSOMhist.col(uhistIDi);
		if (uSOM_Vti == Eigen::VectorXd::Zero(6)) continue;
		
		// Get linear control actions (displacements)
		for (int ucoordi=0; ucoordi<6; ucoordi++) {
			ulin_Vti(ucoordi) = uSOM_Vti(ucoordi) - VSOMstate(SOM_coord_ctrl(ucoordi));
		}
		ulin2_Vti.col(0) << ulin_Vti(0), ulin_Vti(2), ulin_Vti(4);
		ulin2_Vti.col(1) << ulin_Vti(1), ulin_Vti(3), ulin_Vti(5);
		
		// Obtain Vcloth base
		Vcloth_x << VSOMstate(SOM_coord_ctrl(1))-VSOMstate(SOM_coord_ctrl(0)),
			          VSOMstate(SOM_coord_ctrl(3))-VSOMstate(SOM_coord_ctrl(2)),
			          VSOMstate(SOM_coord_ctrl(5))-VSOMstate(SOM_coord_ctrl(4));
		Vcloth_y << -Vcloth_x(1), Vcloth_x(0), 0;                 
		Vcloth_z = Vcloth_x.cross(Vcloth_y);
		
		Vcloth_x = Vcloth_x/Vcloth_x.norm();
		Vcloth_y = Vcloth_y/Vcloth_y.norm();
		Vcloth_z = Vcloth_z/Vcloth_z.norm();
		VRcloth << Vcloth_x, Vcloth_y, Vcloth_z;
		
		// Linear SOM uses local base variables (rot)
		pos_ini_VSOM.col(0) = VSOMstate.segment(0, nSOM*nSOM);
		pos_ini_VSOM.col(1) = VSOMstate.segment(1*nSOM*nSOM, nSOM*nSOM);
		pos_ini_VSOM.col(2) = VSOMstate.segment(2*nSOM*nSOM, nSOM*nSOM);
		vel_ini_VSOM.col(0) = VSOMstate.segment(3*nSOM*nSOM, nSOM*nSOM);
		vel_ini_VSOM.col(1) = VSOMstate.segment(4*nSOM*nSOM, nSOM*nSOM);
		vel_ini_VSOM.col(2) = VSOMstate.segment(5*nSOM*nSOM, nSOM*nSOM);
	
		pos_ini_VSOM_rot = (VRcloth.inverse() * pos_ini_VSOM.transpose()).transpose();
		vel_ini_VSOM_rot = (VRcloth.inverse() * vel_ini_VSOM.transpose()).transpose();
		
		VSOMstate_rot << pos_ini_VSOM_rot.col(0),
				             pos_ini_VSOM_rot.col(1),
				             pos_ini_VSOM_rot.col(2),
				             vel_ini_VSOM_rot.col(0),
				             vel_ini_VSOM_rot.col(1),
				             vel_ini_VSOM_rot.col(2);
				             
		// Rotate control actions from history
		urot2_Vti = VRcloth.inverse() * ulin2_Vti;
		urot_Vti << urot2_Vti.row(0).transpose(),
								urot2_Vti.row(1).transpose(),
								urot2_Vti.row(2).transpose();
				             
		// Simulate a step
		VSOMstate_rot = A_SOM*VSOMstate_rot + B_SOM*urot_Vti + Ts*f_SOM; 
	
		// Convert back to global axes
		pos_nxt_VSOM_rot.col(0) = VSOMstate_rot.segment(0, nSOM*nSOM);
		pos_nxt_VSOM_rot.col(1) = VSOMstate_rot.segment(1*nSOM*nSOM, nSOM*nSOM);
		pos_nxt_VSOM_rot.col(2) = VSOMstate_rot.segment(2*nSOM*nSOM, nSOM*nSOM);
		vel_nxt_VSOM_rot.col(0) = VSOMstate_rot.segment(3*nSOM*nSOM, nSOM*nSOM);
		vel_nxt_VSOM_rot.col(1) = VSOMstate_rot.segment(4*nSOM*nSOM, nSOM*nSOM);
		vel_nxt_VSOM_rot.col(2) = VSOMstate_rot.segment(5*nSOM*nSOM, nSOM*nSOM);
	
		pos_nxt_VSOM = (VRcloth * pos_nxt_VSOM_rot.transpose()).transpose();
		vel_nxt_VSOM = (VRcloth * vel_nxt_VSOM_rot.transpose()).transpose();
		VSOMstate << pos_nxt_VSOM.col(0),
			           pos_nxt_VSOM.col(1),
			           pos_nxt_VSOM.col(2),
			           vel_nxt_VSOM.col(0),
			           vel_nxt_VSOM.col(1),
			           vel_nxt_VSOM.col(2);
  }
  
  // Filter outlying and incorrect mesh data
  Eigen::VectorXd dSOMstate = (VSOMstate - state_SOM).cwiseAbs();
	double mean_dpos = dSOMstate.segment(0,3*nSOM*nSOM).sum()/(3*nSOM*nSOM);
	if (mean_dpos > MaxVDiff) {
	  ROS_WARN_STREAM("Ignoring Vision data: Distance larger than maximum.");
	  return;
	}
	
	// Other filter ideas:
	// - Corner distance specifically
	// - Mean with previous 1-3 instants and add displacement (most recent u)
	// - Get center and xyz offset of nodes, avg with prev 1-3 wrt their ctr, add ctr back
  // Done in Vision:
  // - (Exp) Moving Average (capture time towards the past, more delay)
  // - Gaussian filter with previous/next points (delay 1 step)
	                
	// Update SOM states (weighted average)
	state_SOM = state_SOM*(1-Wv) + VSOMstate*Wv;
	
}


void usomhpReceived(const mpc_pkg::HorizonControls& msg) {
	
	// Extract control actions across all horizon
	uHp_SOM.col(0) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                   (((vector<double>) msg.u1Hp).data(), Hc);
  uHp_SOM.col(1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                   (((vector<double>) msg.u2Hp).data(), Hc);
  uHp_SOM.col(2) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                   (((vector<double>) msg.u3Hp).data(), Hc);
  uHp_SOM.col(3) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                   (((vector<double>) msg.u4Hp).data(), Hc);
  uHp_SOM.col(4) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                   (((vector<double>) msg.u5Hp).data(), Hc);
  uHp_SOM.col(5) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                   (((vector<double>) msg.u6Hp).data(), Hc);
	
	// Reset index
	usomhpID = 0;
	
}




// ---------------------------
// -------- MAIN PROG --------
// ---------------------------

int main(int argc, char **argv) {

	// Initialize the ROS system and become a node.
  ros::init(argc, argv, "mpc_cl_node");
  ros::NodeHandle rosnh;
  
	// Define client & server objects to all services
  ros::ServiceClient clt_shutdown = rosnh.serviceClient<std_srvs::Empty>("node_shutdown");
  ros::service::waitForService("node_shutdown");
  
  // Define Publishers
  ros::Publisher pub_inidata = rosnh.advertise<mpc_pkg::OptiDataQDMC>
                               ("mpc_controller/opti_inidata", 1000);
  ros::Publisher pub_usom = rosnh.advertise<mpc_pkg::TwoPoints>
                            ("mpc_controller/u_SOM", 1000);
  ros::Publisher pub_utcp = rosnh.advertise<geometry_msgs::PoseStamped>
                            ("mpc_controller/u_TCP", 1000);
  ros::Publisher pub_uwam = rosnh.advertise<cartesian_msgs::CartesianCommand>
                            ("iri_wam_controller/CartesianControllerNewGoal", 1000);
  
  // Define Subscribers
  ros::Subscriber sub_somstate = rosnh.subscribe("mpc_controller/state_SOM",
                                                 1000, &somstateReceived);
  ros::Subscriber sub_usomhp = rosnh.subscribe("mpc_controller/u_som_hp",
                                                 1000, &usomhpReceived);
  
  // Get parameters from launch
  if(!rosnh.getParam("/qdmc_rt_node/datapath", datapath)) {
  	ROS_ERROR("Need to define the datapath (where reference csv files are) parameter.");
  }
  if(!rosnh.getParam("/qdmc_rt_node/NTraj", NTraj)) {
  	ROS_ERROR("Need to define the NTraj (reference trajectory number) parameter.");
  }
  if(!rosnh.getParam("/qdmc_rt_node/nSOM", nSOM)) {
  	ROS_ERROR("Need to define the nSOM (SOM mesh side size) parameter.");
  }
  if(!rosnh.getParam("/qdmc_rt_node/nCOM", nCOM)) {
  	ROS_ERROR("Need to define the nCOM (COM mesh side size) parameter.");
  }
  if(!rosnh.getParam("/qdmc_rt_node/Hp", Hp)) {
  	ROS_ERROR("Need to define the Hp (prediction horizon) parameter.");
  }
  if(!rosnh.getParam("/qdmc_rt_node/Hc", Hc)) {
    ROS_ERROR("Need to define the Hc (control horizon) parameter.");
  }
    if(!rosnh.getParam("/qdmc_rt_node/Nss", Nss)) {
      ROS_ERROR("Need to define the Nss (model length) parameter.");
    }
  if(!rosnh.getParam("/qdmc_rt_node/Ts", Ts)) {
  	ROS_ERROR("Need to define the Ts (sample time) parameter.");
  }
  if(!rosnh.getParam("/qdmc_rt_node/Wv", Wv)) {
  	ROS_ERROR("Need to define the Wv (weight of the Vision feedback) parameter.");
  }
  
  
  
 	// --------------------------------
	// 1. INITIAL PARAMETERS DEFINITION
	// --------------------------------
	
	// Load trajectories to follow
  Eigen::MatrixXd phi_l_Traj = getCSVcontent(datapath+"ref_"+to_string(NTraj)+"L.csv");
  Eigen::MatrixXd phi_r_Traj = getCSVcontent(datapath+"ref_"+to_string(NTraj)+"R.csv");
	
	// Get implied cloth size and position
	Eigen::MatrixXd dphi_corners1 = phi_r_Traj.row(0) - phi_l_Traj.row(0);
	double lCloth = dphi_corners1.norm();
  Eigen::Vector3d cCloth;
  cCloth = (phi_r_Traj.row(0) + phi_l_Traj.row(0))/2;
  cCloth(2) += lCloth/2;
  double aCloth = atan2(dphi_corners1(1), dphi_corners1(0));
  
  
  // Initialization of SOM variable
  LinMdl SOM;
  
  // Load parameter table and get row for SOM
  // LUT cols: MdlSz, Ts, theta[7]
  Eigen::MatrixXd thetaLUT = getCSVcontent(datapath+"ThetaMdl_LUT.csv");
  bool SOMparams = false;
  for (int LUTi=0; LUTi<thetaLUT.rows(); LUTi++) {
		if (thetaLUT(LUTi, 1) != Ts) continue;
		
		Eigen::MatrixXd theta_i = thetaLUT.block(LUTi,2, 1,7);
		
		// Set model parameters
		if (thetaLUT(LUTi, 0) == nSOM) {
			if (SOMparams) ROS_WARN_STREAM("Several rows match the same SOM parameters!");
			
		  SOM.stiffness << theta_i(0), theta_i(1), theta_i(2);
  		SOM.damping   << theta_i(3), theta_i(4), theta_i(5);
  		SOM.z_sum      = theta_i(6);
  		SOMparams = true;
		}
		
  }
  // If no row matched the settings
  if (!SOMparams) {
  	ROS_WARN_STREAM("No rows match the SOM parameters! Setting everything to default.");
  	
  	nSOM = 4;
		SOM.stiffness << -200.0, -15.0, -200.0;
  	SOM.damping   << -4.0, -2.5, -4.0;
  	SOM.z_sum      =  0.03;
	}
	
	
	// Rest of SOM Definition. Linear model!
  int nxS  = nSOM;
  int nyS  = nSOM;
  SOM.row  = nxS;
  SOM.col  = nyS;
  SOM.mass = 0.1;
  SOM.grav = 9.8;
  SOM.dt   = Ts;
  int SOMlength = nxS*nyS;
  int COMlength = nCOM*nCOM;

	// Important Coordinates (upper/lower corners in x,y,z)
  SOM_nodes_ctrl << nyS*(nxS-1), nyS*nxS-1;
	SOM_coord_ctrl << SOM_nodes_ctrl(0), SOM_nodes_ctrl(1),
										SOM_nodes_ctrl(0)+nxS*nyS, SOM_nodes_ctrl(1)+nxS*nyS,
										SOM_nodes_ctrl(0)+2*nxS*nyS, SOM_nodes_ctrl(1)+2*nxS*nyS;
  SOM.coord_ctrl = SOM_coord_ctrl;
	Eigen::VectorXd SOM_coord_lc(6);
	SOM_coord_lc << 0, nSOM-1, nSOM*nSOM, nSOM*nSOM+nSOM-1, 2*nSOM*nSOM, 2*nSOM*nSOM+nSOM-1;
	SOM.coord_lc = SOM_coord_lc;
  
  
  // SOM Initialization: Mesh on XZ plane
  Eigen::MatrixXd pos(SOMlength,3);
  pos = create_lin_mesh(lCloth, nSOM, cCloth, aCloth);

  
  // Define initial position of the nodes (for ext_force)
  // Second half of the vector is velocities (initial = 0)
  Eigen::VectorXd x_ini_SOM(2*3*nxS*nyS);
  x_ini_SOM.setZero(2*3*nxS*nyS);

  x_ini_SOM.segment(0, SOMlength) = pos.col(0);
  x_ini_SOM.segment(SOMlength, SOMlength) = pos.col(1);
  x_ini_SOM.segment(2*SOMlength, SOMlength) = pos.col(2);
  
  // Reduce initial SOM position to COM size if necessary
  Eigen::VectorXd reduced_pos(3*COMlength);
  Eigen::VectorXd reduced_vel(3*COMlength);
  tie(reduced_pos, reduced_vel) = take_reduced_mesh(x_ini_SOM.segment(0,3*SOMlength),
  																									x_ini_SOM.segment(3*SOMlength,3*SOMlength),
  																									nSOM, nCOM);
  Eigen::VectorXd x_ini_COM(2*3*COMlength);
  x_ini_COM.setZero(2*3*COMlength);
  x_ini_COM.segment(0,3*COMlength) = reduced_pos;
  
  // Rotate initial COM and SOM positions to XZ plane
  Eigen::Matrix3d RCloth_ini;
  RCloth_ini << cos(aCloth), -sin(aCloth), 0,
  				      sin(aCloth),  cos(aCloth), 0,
  				               0,             0, 1;
  
	Eigen::MatrixXd posSOM_XZ(SOMlength,3);
	posSOM_XZ = (RCloth_ini.inverse() * pos.transpose()).transpose();

  
	// Get the linear model for the SOM
 	A_SOM.resize(6*SOMlength, 6*SOMlength);
 	B_SOM.resize(6*SOMlength, 6);
 	f_SOM.resize(6*SOMlength);
 	A_SOM.setZero(6*SOMlength,6*SOMlength);
 	B_SOM.setZero(6*SOMlength,6);
 	f_SOM.setZero(6*SOMlength);
 	
	tie(SOM, A_SOM, B_SOM, f_SOM) = init_linear_model(SOM, posSOM_XZ);
  
  
  // INITIAL INFO
  ROS_INFO_STREAM("\n- Executing Reference Trajectory: "<<NTraj<<
  								"\n- Reference Trajectory has "<<phi_l_Traj.rows() <<" points."<<endl<<
  								"\n- Sample Time (s): \t"<<Ts<<
  								"\n- Prediction Horizon: \t"<<Hp<<
  								"\n- Model sizes: \tnSOM="<<nSOM<<", nCOM="<<nCOM<<endl<<
									"\n- Cloth Side Length (m): \t "<<lCloth<<
                  "\n- Cloth Initial Center: \t"<<cCloth.transpose()<<
                  "\n- Cloth Initial Angle (rad): \t "<<aCloth<<endl);
  

	// ----------------------------------
	// 2. EXECUTION OF THE REAL TIME LOOP
	// ----------------------------------
	
	// Initial controls
	Eigen::VectorXd u_ini(6);
	Eigen::VectorXd u_bef(6);
	Eigen::VectorXd u_SOM(6);
	
	for (int i=0; i<SOM.coord_ctrl.size(); i++) {
		u_ini(i) = x_ini_SOM(SOM.coord_ctrl(i));
	}
	u_bef = u_ini;
	u_SOM = u_ini;
	uHp_SOM.setZero(Hc,6);
	for (int i=0; i<Hc; i++) {
		uHp_SOM.row(i) = u_SOM.transpose();
	}
	
	// Get cloth orientation (rotation matrix)
	Eigen::Vector3d cloth_x(u_SOM(1)-u_SOM(0),
	                        u_SOM(3)-u_SOM(2),
	                        u_SOM(5)-u_SOM(4));
	Eigen::Vector3d cloth_y(-cloth_x(1), cloth_x(0), 0);
	Eigen::Vector3d cloth_z = cloth_x.cross(cloth_y);
	
	cloth_x = cloth_x/cloth_x.norm();
	cloth_y = cloth_y/cloth_y.norm();
	cloth_z = cloth_z/cloth_z.norm();
	
	Eigen::Matrix3d Rcloth;
	Eigen::Matrix3d Rtcp;
	Rcloth << cloth_x, cloth_y,  cloth_z;
	Rtcp   << cloth_y, cloth_x, -cloth_z;
	
	// Auxiliary variables for base changes
	Eigen::VectorXd phi_red(3*COMlength);
  Eigen::VectorXd dphi_red(3*COMlength);
	Eigen::MatrixXd pos_ini_COM(COMlength,3);
	Eigen::MatrixXd vel_ini_COM(COMlength,3);
	Eigen::MatrixXd pos_ini_COM_rot(COMlength,3);
	Eigen::MatrixXd vel_ini_COM_rot(COMlength,3);
	Eigen::VectorXd x_ini_COM_rot(2*3*COMlength);
	Eigen::MatrixXd pos_ini_SOM(SOMlength,3);
	Eigen::MatrixXd vel_ini_SOM(SOMlength,3);
	Eigen::MatrixXd pos_ini_SOM_rot(SOMlength,3);
	Eigen::MatrixXd vel_ini_SOM_rot(SOMlength,3);
	Eigen::VectorXd state_SOM_rot(2*3*SOMlength);
	Eigen::MatrixXd pos_nxt_SOM(SOMlength,3);
	Eigen::MatrixXd vel_nxt_SOM(SOMlength,3);
	Eigen::MatrixXd pos_nxt_SOM_rot(SOMlength,3);
	Eigen::MatrixXd vel_nxt_SOM_rot(SOMlength,3);
  Eigen::MatrixXd Traj_l_Hp_rot(Hp,3);
  Eigen::MatrixXd Traj_r_Hp_rot(Hp,3);
	Eigen::VectorXd u_lin(6);
	Eigen::MatrixXd u_lin2(3,2);
	Eigen::MatrixXd u_rot2(3,2);
	Eigen::VectorXd u_rot(6);
	Eigen::VectorXd pos_upper_corners(6);
	
	// Reference trajectory in horizon
	Eigen::MatrixXd Traj_l_Hp(Hp,3);
	Eigen::MatrixXd Traj_r_Hp(Hp,3);
	
	// Initialize SOM state
	state_SOM.resize(6*SOMlength);
	state_SOM = x_ini_SOM;
	
	// START SIMULATION LOOP
	ros::Duration(1).sleep();
	int tk = 0;
	ros::Rate rate(1/Ts);
  Eigen::VectorXd f_tilde;
  Eigen::MatrixXd f_i_curr = Eigen::VectorXd::Zero(6 * Nss);
  Eigen::MatrixXd d_model_mismatch = Eigen::VectorXd::Zero(6 * Hp);
  
  Eigen::MatrixXd Psi_full = getCSVcontent(datapath+"Psi.csv");
  Eigen::MatrixXd Gamma_full = getCSVcontent(datapath+"Gamma.csv");
  Eigen::MatrixXd S_tilde_ss_full = getCSVcontent(datapath+"S_tilde_ss.csv");
  Eigen::SparseMatrix<double> Psi = Psi_full.sparseView();
  Eigen::SparseMatrix<double> Gamma = Gamma_full.sparseView();
  Eigen::SparseMatrix<double> S_tilde_ss = S_tilde_ss_full.sparseView();
    
	while(rosnh.ok()) {
	
		// Check subscriptions
		ros::spinOnce();
	
		// Slice reference (constant in the window near the end)
	  if(tk >= phi_l_Traj.rows()-(Hp+1)) {
	    Traj_l_Hp = Eigen::VectorXd::Ones(Hp+1)*phi_l_Traj.bottomRows(1);
	    Traj_r_Hp = Eigen::VectorXd::Ones(Hp+1)*phi_r_Traj.bottomRows(1);
	  }
	  else{
	    Traj_l_Hp = phi_l_Traj.block(tk,0, Hp+1,3);
	    Traj_r_Hp = phi_r_Traj.block(tk,0, Hp+1,3);
	  }
	  
	  // Get reduced states (SOM->COM)
    tie(phi_red, dphi_red) = take_reduced_mesh(state_SOM.segment(0,3*SOMlength),
  																						 state_SOM.segment(3*SOMlength,3*SOMlength),
  																						 nSOM, nCOM);
  	
  	// Update COM state																					 
  	x_ini_COM << phi_red, dphi_red;
		// Rotate initial position to cloth base
		pos_ini_COM.col(0) = x_ini_COM.segment(0, COMlength);
		pos_ini_COM.col(1) = x_ini_COM.segment(1*COMlength, COMlength);
		pos_ini_COM.col(2) = x_ini_COM.segment(2*COMlength, COMlength);
		vel_ini_COM.col(0) = x_ini_COM.segment(3*COMlength, COMlength);
		vel_ini_COM.col(1) = x_ini_COM.segment(4*COMlength, COMlength);
		vel_ini_COM.col(2) = x_ini_COM.segment(5*COMlength, COMlength);
		
		pos_ini_COM_rot = (Rcloth.inverse() * pos_ini_COM.transpose()).transpose();
		vel_ini_COM_rot = (Rcloth.inverse() * vel_ini_COM.transpose()).transpose();
		x_ini_COM_rot << pos_ini_COM_rot.col(0),
		                 pos_ini_COM_rot.col(1),
		                 pos_ini_COM_rot.col(2),
		                 vel_ini_COM_rot.col(0),
		                 vel_ini_COM_rot.col(1),
		                 vel_ini_COM_rot.col(2);
    for (int i=0; i<SOM.coord_ctrl.size(); i++) {
		    pos_upper_corners(i) = x_ini_COM_rot(SOM.coord_ctrl(i));
	  }
    if (tk == 0){
      Eigen::VectorXd x_free_response = x_ini_COM_rot;  
      for(int j = 0; j < Nss; ++j){
        for(int i = 0; i < 6; ++i){
            f_i_curr(i * Nss + j) = x_ini_COM_rot(SOM.coord_lc(i));            
        }
        //ROS_WARN_STREAM("-------" << std::endl << f_i_curr << std::endl);
        //sleep(5);  
        x_free_response = A_SOM * x_free_response + Ts * f_SOM;
      }
    }
   
    f_tilde = Psi * f_i_curr + d_model_mismatch;                
		Traj_l_Hp_rot = (Rcloth.inverse() * Traj_l_Hp.transpose()).transpose();
		Traj_r_Hp_rot = (Rcloth.inverse() * Traj_r_Hp.transpose()).transpose();

		
		// Publish x_ini_COM_rot, u_rot, Traj_x_Hp_rot, Rcloth and u_bef to optimizer
		mpc_pkg::OptiDataQDMC optiData_pub;
		optiData_pub.f_tilde = vector<double>(f_tilde.data(),
												 f_tilde.size()+f_tilde.data());
		optiData_pub.u_rot = vector<double>(pos_upper_corners.data(), pos_upper_corners.size()+pos_upper_corners.data());
		optiData_pub.traj_left_x = vector<double>(Traj_l_Hp_rot.col(0).data(),
									Traj_l_Hp_rot.col(0).size()+Traj_l_Hp_rot.col(0).data());
		optiData_pub.traj_left_y = vector<double>(Traj_l_Hp_rot.col(1).data(),
									Traj_l_Hp_rot.col(1).size()+Traj_l_Hp_rot.col(1).data());
		optiData_pub.traj_left_z = vector<double>(Traj_l_Hp_rot.col(2).data(),
									Traj_l_Hp_rot.col(2).size()+Traj_l_Hp_rot.col(2).data());
		optiData_pub.traj_right_x = vector<double>(Traj_r_Hp_rot.col(0).data(),
									 Traj_r_Hp_rot.col(0).size()+Traj_r_Hp_rot.col(0).data());
		optiData_pub.traj_right_y = vector<double>(Traj_r_Hp_rot.col(1).data(),
									 Traj_r_Hp_rot.col(1).size()+Traj_r_Hp_rot.col(1).data());
		optiData_pub.traj_right_z = vector<double>(Traj_r_Hp_rot.col(2).data(),
									 Traj_r_Hp_rot.col(2).size()+Traj_r_Hp_rot.col(2).data());
		optiData_pub.Rcloth_x = vector<double>(cloth_x.data(), cloth_x.size()+cloth_x.data());
		optiData_pub.Rcloth_y = vector<double>(cloth_y.data(), cloth_y.size()+cloth_y.data());
		optiData_pub.Rcloth_z = vector<double>(cloth_z.data(), cloth_z.size()+cloth_z.data());
		optiData_pub.u_bef = vector<double>(u_bef.data(), u_bef.size()+u_bef.data());
		
		pub_inidata.publish(optiData_pub);
		
		
		// Get corresponding u_SOM from array and index
		u_SOM = uHp_SOM.row(usomhpID).transpose();
		u_lin = u_SOM - u_bef;
		usomhpID = min(usomhpID+1, Hc-1);
		

		// Update uSOM history (6xNH)
		uSOMhist.col(uhistID) = u_SOM;
		uhistID = (uhistID+1)%MaxVSteps;
		
		
		// Linear SOM uses local base variables (rot)
		// - Rotate state vector
		pos_ini_SOM.col(0) = state_SOM.segment(0, SOMlength);
		pos_ini_SOM.col(1) = state_SOM.segment(1*SOMlength, SOMlength);
		pos_ini_SOM.col(2) = state_SOM.segment(2*SOMlength, SOMlength);
		vel_ini_SOM.col(0) = state_SOM.segment(3*SOMlength, SOMlength);
		vel_ini_SOM.col(1) = state_SOM.segment(4*SOMlength, SOMlength);
		vel_ini_SOM.col(2) = state_SOM.segment(5*SOMlength, SOMlength);
		
		pos_ini_SOM_rot = (Rcloth.inverse() * pos_ini_SOM.transpose()).transpose();
		vel_ini_SOM_rot = (Rcloth.inverse() * vel_ini_SOM.transpose()).transpose();
		state_SOM_rot << pos_ini_SOM_rot.col(0),
		                 pos_ini_SOM_rot.col(1),
		                 pos_ini_SOM_rot.col(2),
		                 vel_ini_SOM_rot.col(0),
		                 vel_ini_SOM_rot.col(1),
		                 vel_ini_SOM_rot.col(2);
		                 
		// - Rotate control input
		u_lin2.col(0) << u_lin(0), u_lin(2), u_lin(4);
		u_lin2.col(1) << u_lin(1), u_lin(3), u_lin(5);
		u_rot2 = Rcloth.inverse() * u_lin2;
		u_rot << u_rot2.row(0).transpose(), u_rot2.row(1).transpose(), u_rot2.row(2).transpose();
		

		// Simulate a step
	  state_SOM_rot = A_SOM*state_SOM_rot + B_SOM*u_rot + Ts*f_SOM;
	  
	  
	  // Convert back to global axes
	  pos_nxt_SOM_rot.col(0) = state_SOM_rot.segment(0, SOMlength);
		pos_nxt_SOM_rot.col(1) = state_SOM_rot.segment(1*SOMlength, SOMlength);
		pos_nxt_SOM_rot.col(2) = state_SOM_rot.segment(2*SOMlength, SOMlength);
		vel_nxt_SOM_rot.col(0) = state_SOM_rot.segment(3*SOMlength, SOMlength);
		vel_nxt_SOM_rot.col(1) = state_SOM_rot.segment(4*SOMlength, SOMlength);
		vel_nxt_SOM_rot.col(2) = state_SOM_rot.segment(5*SOMlength, SOMlength);
		
		pos_nxt_SOM = (Rcloth * pos_nxt_SOM_rot.transpose()).transpose();
		vel_nxt_SOM = (Rcloth * vel_nxt_SOM_rot.transpose()).transpose();
		state_SOM << pos_nxt_SOM.col(0),
		             pos_nxt_SOM.col(1),
		             pos_nxt_SOM.col(2),
		             vel_nxt_SOM.col(0),
		             vel_nxt_SOM.col(1),
		             vel_nxt_SOM.col(2);
		             
		// Update u_bef
		for (int i=0; i<SOM.coord_ctrl.size(); i++) {
			u_bef(i) = state_SOM(SOM.coord_ctrl(i));
		}
  	
  	// Get new Cloth orientation (rotation matrix)
  	cloth_x << u_SOM(1)-u_SOM(0),
	             u_SOM(3)-u_SOM(2),
	             u_SOM(5)-u_SOM(4);
		cloth_y << -cloth_x(1), cloth_x(0), 0;
		cloth_z = cloth_x.cross(cloth_y);
	
		cloth_x = cloth_x/cloth_x.norm();
		cloth_y = cloth_y/cloth_y.norm();
		cloth_z = cloth_z/cloth_z.norm();
		
		Eigen::MatrixXd old_Rcloth = Rcloth;
	
		Rcloth << cloth_x, cloth_y,  cloth_z;
		Rtcp   << cloth_y, cloth_x, -cloth_z;
  	
  	//ROS_WARN_STREAM(endl<<Rcloth <<endl << old_Rcloth  <<endl);
  	
  	// Transform orientation to quaternion
  	tf2::Matrix3x3 tfRtcp;
		tfRtcp.setValue(Rtcp(0,0), Rtcp(0,1), Rtcp(0,2),
		              	Rtcp(1,0), Rtcp(1,1), Rtcp(1,2),
		              	Rtcp(2,0), Rtcp(2,1), Rtcp(2,2));
		
		tf2::Quaternion tfQtcp;
		tfRtcp.getRotation(tfQtcp);
		geometry_msgs::Quaternion Qtcp = tf2::toMsg(tfQtcp);
		// Publish control action: Two corners
		mpc_pkg::TwoPoints u_SOM_pub;
		
		u_SOM_pub.pt1.x = u_SOM(0);
		u_SOM_pub.pt2.x = u_SOM(1);
		u_SOM_pub.pt1.y = u_SOM(2);
		u_SOM_pub.pt2.y = u_SOM(3);
		u_SOM_pub.pt1.z = u_SOM(4);
		u_SOM_pub.pt2.z = u_SOM(5);
		
		pub_usom.publish(u_SOM_pub);
		
		// Cloth-TCP Base offset in Robot coords
		Eigen::Vector3d tcp_offset = Rcloth * tcp_offset_local;
		
		// Publish control action (TCP) if not NaN
		if (!u_SOM.hasNaN()) {
		  cartesian_msgs::CartesianCommand u_WAM_pub;
		  geometry_msgs::PoseStamped u_TCP_pub;
		
		  u_TCP_pub.header.stamp = ros::Time::now();
		  u_TCP_pub.header.frame_id = "iri_wam_link_base";
		  u_TCP_pub.pose.position.x = (u_SOM(0)+u_SOM(1))/2 + tcp_offset(0);
		  u_TCP_pub.pose.position.y = (u_SOM(2)+u_SOM(3))/2 + tcp_offset(1);
		  u_TCP_pub.pose.position.z = (u_SOM(4)+u_SOM(5))/2 + tcp_offset(2);
		  u_TCP_pub.pose.orientation = Qtcp;
		
		  u_WAM_pub.desired_pose = u_TCP_pub;
		  u_WAM_pub.duration = 0.010;
		
		  pub_utcp.publish(u_TCP_pub);
		  pub_uwam.publish(u_WAM_pub);
		}
	  Eigen::VectorXd u_rot_padded = Eigen::VectorXd::Zero(6 * Hc);
		for(int i = 0; i < 6; ++i){
        u_rot_padded(i * Hc) = u_rot(i);
    }
		f_i_curr = Gamma * f_i_curr + S_tilde_ss * u_rot_padded;
        
    //Update model mismatch and free response array
    Eigen::VectorXd meas_output = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd meas_output_rot = Eigen::VectorXd::Zero(6);
    for(int i = 0; i < 6; ++i){
        meas_output(i) = state_SOM(SOM.coord_lc(i));
        /*if(i == 0){
          ROS_INFO_STREAM(meas_output(i));
        }*/
    } 

    Eigen::VectorXd meas_output_rot_l = (Rcloth.inverse() * meas_output.block(0, 0, 3, 1));
    Eigen::VectorXd meas_output_rot_r = (Rcloth.inverse() * meas_output.block(3, 0, 3, 1));
    meas_output_rot.block(0, 0, 3, 1) = meas_output_rot_l;
    meas_output_rot.block(3, 0, 3, 1) = meas_output_rot_r;
    for(int i = 0; i < 6; ++i){
        double fkk = f_i_curr(i * Nss);
        double model_mismatch = meas_output_rot(i) - fkk;
        /*if(i == 0){
          ROS_WARN_STREAM(meas_output_rot(i) << " " << model_mismatch);
        }*/

        for(int j = 0; j < Hp; ++j){
            d_model_mismatch(i * Hp + j) = model_mismatch;
        }
    }
    //ROS_WARN_STREAM(u_SOM);
    //sleep(3);
    //sleep(2);

		// Increase counter
		tk++;
		
		
		// End after trajectory is completed
		if (tk > phi_l_Traj.rows() + 4/Ts) {
		  break;
		}
		
		// Execute at a fixed rate
		rate.sleep();
	
	}
	// END LOOP
	
	// Shutdown node and optimizer
	std_srvs::Empty::Request call_req;
	std_srvs::Empty::Response call_resp;
	bool shutdown_success = clt_shutdown.call(call_req, call_resp);
  ROS_INFO_STREAM("Reached the end!" << endl);
  
}































