#include "./SINS.h"

int main(){
    SINS::SINS_base sins_base;
   /*  Vector3d a(1,1,1);
    Matrix3d b;
    b = sins_base.askew(a);
    cout<<" b is"<<b<<endl; */
    
    SINS::SINS_KF sins_kf;
    /* sins_kf.readIMUdata("/home/e/图片/vins_data/imu.txt");
    sins_kf.readXYZdata("/home/e/图片/vins_data/XYZ.txt");
     */
    sins_kf.SINS_KF::kfmain();
    //cout<<sins_kf.Gyro[5]<<endl;
    return 0;
}