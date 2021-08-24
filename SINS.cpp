#include "./SINS.h"
namespace SINS{
SINS_base::SINS_base(/* args */)
{
}

SINS_base::~SINS_base()
{
}
//三维向量的反对称矩阵
Matrix3d SINS_base::askew(Vector3d v){
    Matrix3d out;
    out<<0,-v(2),v(1),
        v(2),0,-v(0),
        -v(1),v(0),0;
    return out;
}
//姿态角转化为姿态阵
Matrix3d SINS_base::a2mat(Vector3d v){
    Matrix3d out;
    double si,sj,sk,ci,cj,ck;
    si = sin(v(0));sj = sin(v(1));sk = sin(v(2));
    ci = cos(v(0));cj = cos(v(1));ck = cos(v(2));
    out<< cj*ck-si*sj*sk , -ci*sk , sj*ck+si*cj*sk,
          cj*sk+si*sj*ck , ci*ck , sj*sk-si*cj*ck,
          -ci*sj , si , ci*cj;
    return out;
}
//姿态阵转化为姿态角
Vector3d SINS_base::m2att(Matrix3d m){
    Vector3d out;
    if(abs(m(2,1))<= 0.999999){
        out << asin(m(2,1)) , -atan2(m(2,0),m(2,2)) , -atan2(m(0,1),m(1,1));
    }
    else{
        out << asin(m(2,1)) , atan2(m(0,2),m(0,0)) , 0;
    }
    return out;
}
//姿态角转化为四元数
Quaterniond SINS_base::a2qua(Vector3d v){
    Quaterniond out;
    double si,sj,sk,ci,cj,ck;
    si = sin(v(0)/2);sj = sin(v(1)/2);sk = sin(v(2)/2);
    ci = cos(v(0)/2);cj = cos(v(1)/2);ck = cos(v(2)/2);
    out.w() = ci*cj*ck-si*sj*sk;
    out.x() = si*cj*ck-ci*sj*sk,
    out.y() = ci*sj*ck+si*cj*sk,
    out.z() = ci*cj*sk+si*sj*ck;
    return out;
}
//四元数转化为姿态角
Vector3d SINS_base::q2att(Quaterniond q){
    Vector3d out;
    out = m2att(q2mat(q));
    return out;
}
//姿态阵转化为四元数   
Quaterniond SINS_base::m2qua(Matrix3d m){
    Quaterniond out;
    double q0,q1,q2,q3;
    double C11=m(0,0);double C12=m(0,1);double C13=m(0,2);
    double C21=m(1,0);double C22=m(1,1);double C23=m(1,2);
    double C31=m(2,0);double C32=m(2,1);double C33=m(2,2);
    if(C11>=C22+C33){
        q1 = 0.5*sqrt(1+C11-C22-C33);
        q0 = (C32-C23)/(4*q1); q2 = (C12+C21)/(4*q1); q3 = (C13+C31)/(4*q1);
    }
    else if(C22>=C11+C33){
        q2 = 0.5*sqrt(1-C11+C22-C33);
        q0 = (C13-C31)/(4*q2); q1 = (C12+C21)/(4*q2); q3 = (C23+C32)/(4*q2);
    }
    else if(C33>=C11+C22){
        q3 = 0.5*sqrt(1-C11-C22+C33);
        q0 = (C21-C12)/(4*q3); q1 = (C13+C31)/(4*q3); q2 = (C23+C32)/(4*q3);
    }    
    else{
        q0 = 0.5*sqrt(1+C11+C22+C33);
        q1 = (C32-C23)/(4*q0); q2 = (C13-C31)/(4*q0); q3 = (C21-C12)/(4*q0);
    }
    out.w() = q0;
    out.x() = q1;
    out.y() = q2;
    out.z() = q3;
    return out;
}
//四元数转化为姿态阵    
Matrix3d SINS_base::q2mat(Quaterniond q){
    Matrix3d out;
    double q11=q.w()*q.w(); double q12=q.w()*q.x(); double q13=q.w()*q.y(); double q14=q.w()*q.z();
    double q22=q.x()*q.x(); double q23=q.x()*q.y(); double q24=q.x()*q.z();
    double q33=q.y()*q.y(); double q34=q.y()*q.z();
    double q44=q.z()*q.z();
    out<< q11+q22-q33-q44, 2*(q23-q14), 2*(q24+q13),
          2*(q23+q14), q11-q22+q33-q44, 2*(q34-q12),
          2*(q24-q13), 2*(q34+q12), q11-q22-q33+q44;
    return out;
}
//旋转矢量转化为变换矩阵
Matrix3d SINS_base::rv2m(Vector3d v){
    Matrix3d out;
    double nm2= v.transpose()*v;// squaredNorm();
    double a,b;
    if(nm2 <1e-8){
        a = 1-nm2*(1/6-nm2/120); b = 0.5-nm2*(1/24-nm2/720);
    }
    else{
        double nm = sqrt(nm2);
        a = sin(nm)/nm ; b = (1-cos(nm))/nm2;
    }
    Matrix3d VX =askew(v);
    out = Eigen::Matrix3d::Identity() + a*VX + b*VX*VX;
    return out;
}
//旋转矢量转化为变化四元数    
Quaterniond SINS_base::rv2q(Vector3d v){
    Quaterniond out;
    double nm2= v.transpose()*v;
    double q0,s;
    if(nm2 <1e-8){
        q0 = 1-nm2*(1/8-nm2/384); s = 0.5-nm2*(1/48-nm2/3840);
    }
    else{
        double nm = sqrt(nm2);
        q0 = cos(nm/2) ; s = sin(nm/2)/nm;
    }
    out.w() = q0;
    out.x() = s*v(0);
    out.y() = s*v(1);
    out.z() = s*v(2);
    return out;
}
//变化四元数转化为旋转矢量
Vector3d SINS_base::q2rv(Quaterniond q){
    Vector3d out;
    if(q.w()<0){
        q.w()=-q.w();q.x()=-q.x();q.y()=-q.y();q.z()=-q.z();
    }
        
    double nmhalf=acos(q.w());//等效旋转矢量模值的一半
    double b;
    if(nmhalf>1e-20){
        b = 2*nmhalf/sin(nmhalf);
    }
    else{
        b = 2;
    }
    out(0) = b*q.x();
    out(1) = b*q.y();
    out(2) = b*q.z();
    return out;
}
//四元数共轭    
Quaterniond SINS_base::qconj(Quaterniond q_in){
    Quaterniond out;
    out.w() = q_in.w();
    out.x() = -q_in.x();
    out.y() = -q_in.y();
    out.z() = -q_in.z();
    //out = q_in.conjugate();
    return out;
}
//四元数归一化
Quaterniond SINS_base::qnormlz(Quaterniond q_in){
    Quaterniond out;
    out = q_in.normalized();
    return out;
}
//四元数相乘    
Quaterniond SINS_base::qmul(Quaterniond q_1,Quaterniond q_2){
    Quaterniond out(
    q_1.w() * q_2.w() - q_1.x() * q_2.x() - q_1.y() * q_2.y() - q_1.z() * q_2.z(),
    q_1.w() * q_2.x() + q_1.x() * q_2.w() + q_1.y() * q_2.z() - q_1.z() * q_2.y(),
    q_1.w() * q_2.y() + q_1.y() * q_2.w() + q_1.z() * q_2.x() - q_1.x() * q_2.z(),
    q_1.w() * q_2.z() + q_1.z() * q_2.w() + q_1.x() * q_2.y() - q_1.y() * q_2.x()
    );
    return out;
}
//四元数乘向量    
Vector3d SINS_base::qmulv(Quaterniond q_1,Vector3d v_1){
    Quaterniond out_q;
    Quaterniond tmp_q(0,v_1(0),v_1(1),v_1(2));
    out_q = qmul(qmul(q_1,tmp_q),qconj(q_1));
    Vector3d out(out_q.x(),out_q.y(),out_q.z());
    return out;
}
//四元数加失准角误差
Quaterniond SINS_base::qaddphi(Quaterniond q,Vector3d phi){
    Quaterniond out;
    out = qmul(rv2q(-phi),q);
    return out;
}
//四元数减失准角误差
Quaterniond SINS_base::qdelphi(Quaterniond q,Vector3d phi){
    Quaterniond out;
    out = qmul(rv2q(phi),q);
    return out;
}
//由计算四元数和真实四元数计算失准角误差  
Vector3d SINS_base::qq2phi(Quaterniond q_1,Quaterniond q_2){
    Quaterniond q_err;
    q_err =qmul(q_2,qconj(q_1));
    Vector3d out;
    out = q2rv(q_err);
    return out;
}

//圆锥误差划船误差补偿
void SINS_base::cnscl(vector<Vector3d> wm,vector<Vector3d> vm,Vector3d& phim,Vector3d& dvbm){
    /* 分别定义两种情况 2子样，3子样 */
    Vector3d wmm = Vector3d::Zero();
    Vector3d vmm = Vector3d::Zero();
    Vector3d dphim = Vector3d::Zero();
    Vector3d scullm = Vector3d::Zero();
    //if(wm.size()==2){    
        //for(int i=0;i<wm.size();i++){
        wmm=wm[0]+wm[1];
        vmm=vm[0]+vm[1];
        //}
        Vector3d csw= 2/3*wm[0];
        Vector3d csv= 2/3*vm[0];
        dphim = askew(csw)*wm[1];//圆锥补偿量
        scullm = askew(csw)*vm[1] + askew(csv)*wm[1];//圆锥补偿量
        phim = wmm+dphim;
        dvbm = vmm+0.5*(askew(wmm)*vmm)+scullm;
    //}
   /*  else if(wm.size()==3){
        //暂时用不到先不定义
    }
    else{
        cout <<" 大于3子样，请给合适大小！"<<endl;
    } */
}

/* 该版本为简易版的惯性更新算法 */
void SINS_base::insupdate(Quaterniond& qnb,Vector3d& vn,Vector3d& pos,
            vector<Vector3d> wm,vector<Vector3d> vm,double ts){
    double n =2.0;
    double nts = n*ts;
    Vector3d phim = Vector3d::Zero();
    Vector3d dvbm=Vector3d::Zero();
    cnscl(wm,vm,phim,dvbm);//圆锥误差/划船误差补偿
    Vector3d gcc(0 ,0 ,-9.780325);//固定的重力
    Vector3d vn1 = Vector3d::Zero();
    /* cout<<"phim "<<phim.transpose()<<endl;
    cout<<"dvbm "<<dvbm.transpose()<<endl; */
    vn1 = vn + qmulv(qnb,dvbm) + gcc*nts;  // 速度更新
    vn = (vn+vn1)/2;
    pos = pos + vn*nts;  //使用平均速度对位置进行基于xyz的更新，4.1.60
    vn = vn1;
    qnb = qmul(qnb, rv2q(phim));//姿态更新4.1.8
    qnb = qnormlz(qnb);
    /* cout<<wm[0].transpose()<<" "<<vn1.transpose()<<endl; */
}


////////////////////////////////////////////
/* SINS_earth::SINS_earth()
{
}

SINS_earth::~SINS_earth()
{
}

*/
////////////////////////////////////////////

SINS_KF::SINS_KF()
{
}

SINS_KF::~SINS_KF()
{
} 

//15维状态误差转移矩阵
 Matrix15d SINS_KF::kfft15(Matrix3d Cnb,Vector3d fb ,SINS_base sins_base){
    Matrix3d Mva,Mpv,O33;
    //SINS_base sins_base;
    Mva = sins_base.askew(Cnb*fb);
    Mpv = Matrix3d::Identity();
    //O33 = Matrix3d::Zero();
    Matrix15d Ft = Matrix15d::Zero();
    Ft.block<3,3>(0,9) = -Cnb;
    Ft.block<3,3>(3,0) = Mva;
    Ft.block<3,3>(3,12) = Cnb;
    Ft.block<3,3>(6,3) = Mpv;
    /* Ft = 
    [ O33    O33    O33    -Cnb     O33 
    Mva    O33    O33     O33     Cnb 
    O33    Mpv    O33     O33     O33
    zeros(6,15) ]; cout<<Ft<<endl<<endl;*/
    return Ft;
}

void SINS_KF::kfinit(Matrix15d Qk,Matrix6d Rk,Matrix15d P0, Matrix15d Phikk_1, Matrix6_15d Hk,Matrix15d Gammak){
    /* const int m,n;
    m = Hk.rows(); n = Hk.cols(); */
    //初始化的参数
    kf_Qk = Qk;
    kf_Rk = Rk;
    kf_Pk = P0;
    kf_Xk = Eigen::Matrix<double,15,1>::Zero();
    kf_Phikk_1 = Phikk_1;
    kf_Hk = Hk;
    kf_Gammak = Gammak;
}  
//KALMAN滤波更新
/* flag == 0代表仅时间更新，flag == 1代表量测更新 */
 void SINS_KF::kfupdate(Vector6d Zk,bool flag){
    if(flag == false){
        //cout<<"时间更新!"<<endl;
        kf_Xkk_1 = kf_Phikk_1*kf_Xk;
        kf_Pkk_1 = kf_Phikk_1*kf_Pk*kf_Phikk_1.transpose() + kf_Gammak*kf_Qk*kf_Gammak.transpose();
        kf_Xk = kf_Xkk_1;
        kf_Pk = kf_Pkk_1;
    }
    else if(flag == true){
        //cout<<"量测更新!"<<endl;
        //cout<<"矫正前 Xk"<<kf_Xk.block<6,1>(3,0).transpose()<<endl;
        /* kf_Xkk_1 = kf_Xk;
        kf_Pkk_1 = kf_Pk; */
        /* kf_PXZkk_1 = kf_Pkk_1*kf_Hk.transpose();
        kf_PZkk_1 = kf_Hk*kf_PXZkk_1+kf_Rk;
        kf_Kk = kf_PXZkk_1*kf_PZkk_1.inverse(); */
        kf_Kk = kf_Pkk_1*kf_Hk.transpose()*(kf_Hk*kf_Pkk_1*kf_Hk.transpose()+kf_Rk).inverse();
        kf_Xk = kf_Xkk_1 + kf_Kk*(Zk-kf_Hk*kf_Xkk_1);
        //cout<<"矫正 Xk"<<kf_Xk.block<6,1>(3,0).transpose()<<endl;
        //kf_Pk = kf_Pkk_1- kf_Kk*kf_PZkk_1*kf_Kk.transpose();
        kf_Pk = kf_Pkk_1- kf_Kk*(kf_Hk*kf_Pkk_1*kf_Hk.transpose()+kf_Rk)*kf_Kk.transpose();

    }
    kf_Pk = (kf_Pk+kf_Pk.transpose())/2;//p阵对称化
} 

void SINS_KF::kfmain(){
double nn = 2.0; 
double ts = 1.0/200.0;  
double nts = nn*ts; //子样数和采样时间
//初始位置设置
SINS_base sins_base;
Quaterniond qnb0(0.238261,-0.75761,-0.348629,-0.497711);
Vector3d att0 = sins_base.q2att(qnb0);
Vector3d vn0(0.006246,-0.001431,-0.004391); 
Vector3d pos0(4.460675,-1.680515,0.579614);
Quaterniond qnb;Vector3d vn;Vector3d pos;
qnb = qnb0; vn = vn0; pos = pos0; //姿态、速度和位置初始化
/* 载入数据 读入txt*/
readIMUdata("/home/e/图片/vins_data/imu.txt");
readXYZdata("/home/e/图片/vins_data/XYZ.txt");
//TO DO
//////////////////////////////////////
//添加失准角误差
//Vector3d phi(5.0,5.0,5.0);  qnb = sins_base.qaddphi(qnb, phi);  
Vector3d web(-0.001806,0.02094,0.07687); //设定陀螺常值零偏
Vector3d wdb (-0.020544,0.124837,0.0618);  //加速度计常值零偏
//设置噪声得协方差阵R,Q
Matrix15d Qk = Matrix15d::Zero();
Matrix3d tmp ;
tmp<<web(0)*web(0),0,0,
    0,web(1)*web(1),0,
    0,0,web(2)*web(2);
Qk.block<3,3>(0,0) = tmp*nts;
tmp<<wdb(0)*wdb(0),0,0,
    0,wdb(1)*wdb(1),0,
    0,0,wdb(2)*wdb(2);
Qk.block<3,3>(3,3) = tmp*nts;

Eigen::Matrix<double,6,1> rk;
rk<<0.1,0.1,0.1,                        
    1.0,1.0,1.0; 
Matrix6d Rk = Matrix6d::Zero();
Rk(0,0)=rk(0)*rk(0);Rk(1,1)=rk(1)*rk(1);Rk(2,2)=rk(2)*rk(2);
Rk(3,3)=rk(3)*rk(3);Rk(4,4)=rk(4)*rk(4);Rk(5,5)=rk(5)*rk(5);

//设置量测矩阵
Matrix6_15d Hk = Matrix6_15d::Zero();
Hk.block<6,6>(0,3) = Matrix6d::Identity();
//设置初始协防差
Matrix15d P0 = Matrix15d::Identity();
//kf滤波器初始化
kfinit(Qk, Rk, P0, Matrix15d::Zero(), Hk, Matrix15d::Identity());
//设置仿真时间 
int len = Acc.size();//总步长

/* vector<Vector3d> pos_ans;
vector<Vector6d> bias;
vector<Matrix15d> P_ans;
 */
vector<Vector3d> wm1,vm1;
//nts=0.01s,等价与kk=1,如果
double t = 0.0; //记录导航结果
int kk = 0;
const double update_period = 0.02;
const int update_times = update_period/nts;
for(int k=0;k<10000;k=k+2){
    t+=nts;
    kk+=1;
    wm1.clear();
    vm1.clear();
    //二子样，所以每次载入两个数据
    wm1.push_back(Gyro[k]);wm1.push_back(Gyro[k+1]);
    vm1.push_back(Acc[k]);vm1.push_back(Acc[k+1]);
    //cout<<vm1[0].transpose()<<" , "<<vm1[1].transpose()<<endl;
    sins_base.insupdate(qnb, vn, pos, wm1, vm1, ts);
    //kalman 时间更新&量测更新
    kf_Phikk_1 = Matrix15d::Identity() + 
                kfft15(sins_base.q2mat(qnb),(vm1[0]+vm1[1])/nts,sins_base)*nts;
    kfupdate(Vector6d::Zero(),false);//此时只做时间更新
    
    if (fmod(kk,update_times)<1e-3){
        if(8000>k||k>9000){
        Eigen::Matrix<double,6,1> measure_v_p;
        measure_v_p = XYZ[k];//+ rk.*randn(6,1);
        Eigen::Matrix<double,6,1> estimate_v_p;
        //estimate_v_p = kf_Xk.block<6,1>(3,0);
        estimate_v_p.block<3,1>(0,0) = vn; 
        estimate_v_p.block<3,1>(3,0) =pos;
        cout<<kk<<"  ";
        //cout<<(estimate_v_p).transpose()<<endl;
        cout<<(estimate_v_p).transpose()-(measure_v_p).transpose()<<endl;
        kfupdate((estimate_v_p-measure_v_p),true);//量测更新
    }
    else{
        cout<<"无量测航迹推算!"<<endl;
        Eigen::Matrix<double,6,1> measure_v_p;
        measure_v_p = XYZ[k];//+ rk.*randn(6,1);
        Eigen::Matrix<double,6,1> estimate_v_p;
        //estimate_v_p = kf_Xk.block<6,1>(3,0);
        estimate_v_p.block<3,1>(0,0) = vn; 
        estimate_v_p.block<3,1>(3,0) =pos;
        cout<<kk<<"  ";
        //cout<<(estimate_v_p).transpose()<<endl;
        cout<<(estimate_v_p).transpose()-(measure_v_p).transpose()<<endl;
        cout<<"----------------------------!"<<endl;
    }
    }
    //反馈
    qnb = sins_base.qdelphi(qnb,kf_Xk.block<3,1>(0,0));    
    vn = vn - kf_Xk.block<3,1>(3,0);
    pos = pos - kf_Xk.block<3,1>(6,0);
    //记录输出
   /*  pos_ans.push_back(pos);
    bias.push_back(kf_Xk.block<6,1>(9,0));
    P_ans.push_back(kf_Pk); */
    //kf_Xk = Vector15d::Zero();//反馈置零
    kf_Xk.block<3,1>(0,0) = Vector3d::Zero();
    kf_Xk.block<3,1>(3,0) = Vector3d::Zero();
    kf_Xk.block<3,1>(6,0) = Vector3d::Zero();

     // write result to file
    ofstream foutC("./pos_ans.txt", ios::app);
    foutC.setf(ios::fixed, ios::floatfield);
    foutC.precision(5);
    foutC <<pos<< endl;
    foutC.close();
    
    ofstream fout2("./bias.txt", ios::app);
    fout2.setf(ios::fixed, ios::floatfield);
    fout2.precision(5);
    fout2 <<kf_Xk.block<6,1>(9,0).transpose()<< endl;
    fout2.close();
} 
    cout<<"函数运行结束"<<endl;
} //end kfmain 

void SINS_KF::readIMUdata(string sIMU_file_dir){
    ifstream fsIMU;
	fsIMU.open(sIMU_file_dir.c_str());
	if (!fsIMU.is_open()){
		cerr << "Failed to open IMU data file! " << sIMU_file_dir << endl;
		return;
	}
    std::string sIMU_line;
	double dStampNSec = 0.0;
    Vector3d vGyro;
	Vector3d vAcc;
	while (std::getline(fsIMU, sIMU_line) && !sIMU_line.empty()){
		std::istringstream ssImuData(sIMU_line); 
        /* 读入数据，单位rad/s , m/s^2 */
 		ssImuData >> dStampNSec >> vGyro.x() >> vGyro.y() >> vGyro.z() >> vAcc.x() >> vAcc.y() >> vAcc.z();
        Gyro.push_back(vGyro);
        Acc.push_back(vAcc);
		usleep(50);
	}
    cout << "IMU数据读取完毕！" << endl;
	fsIMU.close();
}
void SINS_KF::readXYZdata(string sXYZ_file_dir){
    ifstream fsXYZ;
	fsXYZ.open(sXYZ_file_dir.c_str());
	if (!fsXYZ.is_open()){
		cerr << "Failed to open XYZ data file! " << sXYZ_file_dir << endl;
		return;
	}
    std::string sXYZ_line;
	double dStampNSec = 0.0;
    Vector3d vXYZ_vol;
    Vector3d vXYZ_pos;
    Vector6d vXYZ;

	while (std::getline(fsXYZ, sXYZ_line) && !sXYZ_line.empty()){
		std::istringstream ssXYZData(sXYZ_line); 
        /* 读入数据，单位m */
		ssXYZData >> dStampNSec >> vXYZ_vol.x() >> vXYZ_vol.y() >> vXYZ_vol.z()
                        >> vXYZ_pos.x() >> vXYZ_pos.y() >> vXYZ_pos.z();
        vXYZ<<vXYZ_vol.x(),vXYZ_vol.y(),vXYZ_vol.z(),vXYZ_pos.x(),vXYZ_pos.y(),vXYZ_pos.z();
        XYZ.push_back(vXYZ);
		usleep(50);
	}
    cout << "真实轨迹数据读取完毕！" << endl;
	fsXYZ.close();
}  

}//end namespace
