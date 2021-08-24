#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>
#include <fstream>
using namespace Eigen;
using namespace std;
typedef Quaternion<double>  Quaterniond;
typedef Eigen::Matrix<double,15,15> Matrix15d;
typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Matrix<double,6,15> Matrix6_15d;
typedef Eigen::Matrix<double,15,6> Matrix15_6d;
typedef Eigen::Matrix<double,15,1> Vector15d;
typedef Eigen::Matrix<double,6,1> Vector6d;



namespace SINS{
class SINS_base
{
private:
    /* data */
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SINS_base(/* args */);
    ~SINS_base();
    Matrix3d askew(Vector3d v);//三维向量的反对称矩阵

    Matrix3d a2mat(Vector3d v);//姿态角转化为姿态阵
    Vector3d m2att(Matrix3d m);//姿态阵转化为姿态角

    Quaterniond a2qua(Vector3d v);//姿态角转化为四元数
    Vector3d q2att(Quaterniond q);//四元数转化为姿态角

    Quaterniond m2qua(Matrix3d m);//姿态阵转化为四元数
    Matrix3d q2mat(Quaterniond q);//四元数转化为姿态阵

    Matrix3d rv2m(Vector3d v);//旋转矢量转化为变换矩阵

    Quaterniond rv2q(Vector3d v);//旋转矢量转化为变化四元数
    Vector3d q2rv(Quaterniond q);//变化四元数转化为旋转矢量

    Quaterniond qconj(Quaterniond q_in);//四元数共轭
    Quaterniond qnormlz(Quaterniond q_in);//四元数归一化
    Quaterniond qmul(Quaterniond q_1,Quaterniond q_2);//四元数相乘
    Vector3d qmulv(Quaterniond q_1,Vector3d v_1);//四元数乘向量
    Quaterniond qaddphi(Quaterniond q,Vector3d phi);//四元数加失准角误差
    Quaterniond qdelphi(Quaterniond q,Vector3d phi);//四元数减失准角误差
    Vector3d qq2phi(Quaterniond q_1,Quaterniond q_2);//由计算四元数和真实四元数计算失准角误差

    void cnscl(vector<Vector3d> wm,vector<Vector3d> vm,Vector3d& phim,Vector3d& dvbm);//圆锥误差划船误差补偿
    void imuadderr();//TO DO 目前使用的是真实数据暂时不需要添加噪声
    void insupdate(Quaterniond& qnb,Vector3d& vn,Vector3d& pos,vector<Vector3d> wm,vector<Vector3d> vm,double ts);
};

//地球信息，在使用经纬高的时候再更新
/* class SINS_earth
{
private:
    
public:
    SINS_earth();
    ~SINS_earth();
};
*/
class SINS_KF
{
private:
    
public:
    //struct kf_param{
        Matrix15d kf_Qk;//系统噪声的方差
        Matrix6d kf_Rk;//量测噪声的方差
        Matrix15d kf_Pk;//均方误差帧
        Vector15d kf_Xk;//系统状态
        Matrix15d kf_Phikk_1;//状态转移矩阵
        Matrix6_15d kf_Hk;//量测矩阵
        Matrix15d kf_Gammak;
        //
        Vector15d kf_Xkk_1;
        Matrix15d kf_Pkk_1;
        Matrix15_6d kf_PXZkk_1;
        Matrix6d kf_PZkk_1;
        Matrix15_6d kf_Kk;

        vector<Vector3d>  Acc,Gyro;//传感器数据
        vector<Vector6d>  XYZ;//真实轨迹数据

    //};

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SINS_KF();
    ~SINS_KF();
    Matrix15d kfft15(Matrix3d Cnb,Vector3d fb,SINS_base);//15维状态误差转移矩阵
    void kfinit(Matrix15d Qk,Matrix6d Rk,Matrix15d P0, Matrix15d Phikk_1, Matrix6_15d Hk,Matrix15d Gammak);
    void kfupdate(Vector6d Zk,bool flag);
    void kfmain();
    void readIMUdata(string);
    void readXYZdata(string);

};
 



}//end namespace

 
 
 

