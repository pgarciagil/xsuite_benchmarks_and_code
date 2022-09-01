#ifndef XOBJ_STDINT
typedef long int64_t;
typedef char int8_t;
typedef unsigned int uint32_t;
#endif
#ifndef XOBJ_TYPEDEF_ArrNFloat64
#define XOBJ_TYPEDEF_ArrNFloat64
typedef  __global  struct ArrNFloat64_s * ArrNFloat64;
  ArrNFloat64 ArrNFloat64_getp(ArrNFloat64 obj){
  int64_t offset=0;
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ArrNFloat64_len(ArrNFloat64 obj){
  int64_t offset=0;
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ArrNFloat64_get(const ArrNFloat64 obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ArrNFloat64_set(ArrNFloat64 obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ArrNFloat64_getp1(ArrNFloat64 obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_ArrNInt64
#define XOBJ_TYPEDEF_ArrNInt64
typedef  __global  struct ArrNInt64_s * ArrNInt64;
  ArrNInt64 ArrNInt64_getp(ArrNInt64 obj){
  int64_t offset=0;
  return (ArrNInt64)(( __global char*) obj+offset);
}
  int64_t ArrNInt64_len(ArrNInt64 obj){
  int64_t offset=0;
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  int64_t ArrNInt64_get(const ArrNInt64 obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ArrNInt64_set(ArrNInt64 obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=16+i0*8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ArrNInt64_getp1(ArrNInt64 obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_ArrNUint32
#define XOBJ_TYPEDEF_ArrNUint32
typedef  __global  struct ArrNUint32_s * ArrNUint32;
  ArrNUint32 ArrNUint32_getp(ArrNUint32 obj){
  int64_t offset=0;
  return (ArrNUint32)(( __global char*) obj+offset);
}
  int64_t ArrNUint32_len(ArrNUint32 obj){
  int64_t offset=0;
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  uint32_t ArrNUint32_get(const ArrNUint32 obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*4;
  return *( __global uint32_t*)(( __global char*) obj+offset);
}
  void ArrNUint32_set(ArrNUint32 obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=16+i0*4;
  *( __global uint32_t*)(( __global char*) obj+offset)=value;
}
   __global uint32_t* ArrNUint32_getp1(ArrNUint32 obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*4;
  return ( __global uint32_t*)(( __global char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_TriLinearInterpolatedFieldMapData
#define XOBJ_TYPEDEF_TriLinearInterpolatedFieldMapData
typedef  __global  struct TriLinearInterpolatedFieldMapData_s * TriLinearInterpolatedFieldMapData;
  TriLinearInterpolatedFieldMapData TriLinearInterpolatedFieldMapData_getp(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  return (TriLinearInterpolatedFieldMapData)(( __global char*) obj+offset);
}
  double TriLinearInterpolatedFieldMapData_get_x_min(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_x_min(TriLinearInterpolatedFieldMapData obj, double value){
  int64_t offset=0;
  offset+=8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp_x_min(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=8;
  return ( __global double*)(( __global char*) obj+offset);
}
  double TriLinearInterpolatedFieldMapData_get_y_min(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=16;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_y_min(TriLinearInterpolatedFieldMapData obj, double value){
  int64_t offset=0;
  offset+=16;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp_y_min(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=16;
  return ( __global double*)(( __global char*) obj+offset);
}
  double TriLinearInterpolatedFieldMapData_get_z_min(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=24;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_z_min(TriLinearInterpolatedFieldMapData obj, double value){
  int64_t offset=0;
  offset+=24;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp_z_min(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=24;
  return ( __global double*)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_get_nx(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=32;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_nx(TriLinearInterpolatedFieldMapData obj, int64_t value){
  int64_t offset=0;
  offset+=32;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* TriLinearInterpolatedFieldMapData_getp_nx(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=32;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_get_ny(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=40;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_ny(TriLinearInterpolatedFieldMapData obj, int64_t value){
  int64_t offset=0;
  offset+=40;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* TriLinearInterpolatedFieldMapData_getp_ny(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=40;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_get_nz(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=48;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_nz(TriLinearInterpolatedFieldMapData obj, int64_t value){
  int64_t offset=0;
  offset+=48;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* TriLinearInterpolatedFieldMapData_getp_nz(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=48;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  double TriLinearInterpolatedFieldMapData_get_dx(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=56;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_dx(TriLinearInterpolatedFieldMapData obj, double value){
  int64_t offset=0;
  offset+=56;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp_dx(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=56;
  return ( __global double*)(( __global char*) obj+offset);
}
  double TriLinearInterpolatedFieldMapData_get_dy(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=64;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_dy(TriLinearInterpolatedFieldMapData obj, double value){
  int64_t offset=0;
  offset+=64;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp_dy(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=64;
  return ( __global double*)(( __global char*) obj+offset);
}
  double TriLinearInterpolatedFieldMapData_get_dz(const TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=72;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_dz(TriLinearInterpolatedFieldMapData obj, double value){
  int64_t offset=0;
  offset+=72;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp_dz(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=72;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 TriLinearInterpolatedFieldMapData_getp_rho(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=112;
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_len_rho(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=112;
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double TriLinearInterpolatedFieldMapData_get_rho(const TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=112;
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_rho(TriLinearInterpolatedFieldMapData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=112;
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp1_rho(TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=112;
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 TriLinearInterpolatedFieldMapData_getp_phi(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_len_phi(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double TriLinearInterpolatedFieldMapData_get_phi(const TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_phi(TriLinearInterpolatedFieldMapData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp1_phi(TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 TriLinearInterpolatedFieldMapData_getp_dphi_dx(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_len_dphi_dx(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double TriLinearInterpolatedFieldMapData_get_dphi_dx(const TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_dphi_dx(TriLinearInterpolatedFieldMapData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp1_dphi_dx(TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 TriLinearInterpolatedFieldMapData_getp_dphi_dy(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_len_dphi_dy(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double TriLinearInterpolatedFieldMapData_get_dphi_dy(const TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_dphi_dy(TriLinearInterpolatedFieldMapData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp1_dphi_dy(TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 TriLinearInterpolatedFieldMapData_getp_dphi_dz(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t TriLinearInterpolatedFieldMapData_len_dphi_dz(TriLinearInterpolatedFieldMapData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double TriLinearInterpolatedFieldMapData_get_dphi_dz(const TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void TriLinearInterpolatedFieldMapData_set_dphi_dz(TriLinearInterpolatedFieldMapData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* TriLinearInterpolatedFieldMapData_getp1_dphi_dz(TriLinearInterpolatedFieldMapData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_ParticlesData
#define XOBJ_TYPEDEF_ParticlesData
typedef  __global  struct ParticlesData_s * ParticlesData;
  ParticlesData ParticlesData_getp(ParticlesData obj){
  int64_t offset=0;
  return (ParticlesData)(( __global char*) obj+offset);
}
  int64_t ParticlesData_get__capacity(const ParticlesData obj){
  int64_t offset=0;
  offset+=8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__capacity(ParticlesData obj, int64_t value){
  int64_t offset=0;
  offset+=8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp__capacity(ParticlesData obj){
  int64_t offset=0;
  offset+=8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  int64_t ParticlesData_get__num_active_particles(const ParticlesData obj){
  int64_t offset=0;
  offset+=16;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__num_active_particles(ParticlesData obj, int64_t value){
  int64_t offset=0;
  offset+=16;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp__num_active_particles(ParticlesData obj){
  int64_t offset=0;
  offset+=16;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  int64_t ParticlesData_get__num_lost_particles(const ParticlesData obj){
  int64_t offset=0;
  offset+=24;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__num_lost_particles(ParticlesData obj, int64_t value){
  int64_t offset=0;
  offset+=24;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp__num_lost_particles(ParticlesData obj){
  int64_t offset=0;
  offset+=24;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  int64_t ParticlesData_get_start_tracking_at_element(const ParticlesData obj){
  int64_t offset=0;
  offset+=32;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set_start_tracking_at_element(ParticlesData obj, int64_t value){
  int64_t offset=0;
  offset+=32;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp_start_tracking_at_element(ParticlesData obj){
  int64_t offset=0;
  offset+=32;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  double ParticlesData_get_q0(const ParticlesData obj){
  int64_t offset=0;
  offset+=40;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_q0(ParticlesData obj, double value){
  int64_t offset=0;
  offset+=40;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp_q0(ParticlesData obj){
  int64_t offset=0;
  offset+=40;
  return ( __global double*)(( __global char*) obj+offset);
}
  double ParticlesData_get_mass0(const ParticlesData obj){
  int64_t offset=0;
  offset+=48;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_mass0(ParticlesData obj, double value){
  int64_t offset=0;
  offset+=48;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp_mass0(ParticlesData obj){
  int64_t offset=0;
  offset+=48;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_p0c(ParticlesData obj){
  int64_t offset=0;
  offset+=248;
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_p0c(ParticlesData obj){
  int64_t offset=0;
  offset+=248;
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_p0c(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=248;
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_p0c(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=248;
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_p0c(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=248;
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_gamma0(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+56);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_gamma0(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+56);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_gamma0(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+56);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_gamma0(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+56);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_gamma0(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+56);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_beta0(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+64);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_beta0(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+64);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_beta0(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+64);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_beta0(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+64);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_beta0(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+64);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_s(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+72);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_s(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+72);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_s(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+72);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_s(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+72);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_s(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+72);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_x(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_x(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_x(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_x(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_x(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+80);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_y(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_y(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_y(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_y(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_y(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+88);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_px(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_px(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_px(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_px(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_px(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+96);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_py(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_py(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_py(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_py(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_py(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+104);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_zeta(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+112);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_zeta(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+112);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_zeta(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+112);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_zeta(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+112);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_zeta(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+112);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_ptau(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+120);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_ptau(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+120);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_ptau(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+120);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_ptau(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+120);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_ptau(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+120);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_delta(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+128);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_delta(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+128);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_delta(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+128);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_delta(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+128);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_delta(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+128);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_rpp(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+136);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_rpp(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+136);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_rpp(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+136);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_rpp(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+136);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_rpp(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+136);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_rvv(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+144);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_rvv(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+144);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_rvv(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+144);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_rvv(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+144);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_rvv(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+144);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_chi(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+152);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_chi(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+152);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_chi(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+152);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_chi(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+152);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_chi(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+152);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_charge_ratio(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+160);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_charge_ratio(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+160);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_charge_ratio(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+160);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_charge_ratio(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+160);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_charge_ratio(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+160);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNFloat64 ParticlesData_getp_weight(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+168);
  return (ArrNFloat64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_weight(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+168);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  double ParticlesData_get_weight(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+168);
  offset+=16+i0*8;
  return *( __global double*)(( __global char*) obj+offset);
}
  void ParticlesData_set_weight(ParticlesData obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+168);
  offset+=16+i0*8;
  *( __global double*)(( __global char*) obj+offset)=value;
}
   __global double* ParticlesData_getp1_weight(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+168);
  offset+=16+i0*8;
  return ( __global double*)(( __global char*) obj+offset);
}
  ArrNInt64 ParticlesData_getp_particle_id(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+176);
  return (ArrNInt64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_particle_id(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+176);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  int64_t ParticlesData_get_particle_id(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+176);
  offset+=16+i0*8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set_particle_id(ParticlesData obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+176);
  offset+=16+i0*8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp1_particle_id(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+176);
  offset+=16+i0*8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  ArrNInt64 ParticlesData_getp_at_element(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+184);
  return (ArrNInt64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_at_element(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+184);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  int64_t ParticlesData_get_at_element(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+184);
  offset+=16+i0*8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set_at_element(ParticlesData obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+184);
  offset+=16+i0*8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp1_at_element(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+184);
  offset+=16+i0*8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  ArrNInt64 ParticlesData_getp_at_turn(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+192);
  return (ArrNInt64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_at_turn(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+192);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  int64_t ParticlesData_get_at_turn(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+192);
  offset+=16+i0*8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set_at_turn(ParticlesData obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+192);
  offset+=16+i0*8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp1_at_turn(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+192);
  offset+=16+i0*8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  ArrNInt64 ParticlesData_getp_state(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+200);
  return (ArrNInt64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_state(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+200);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  int64_t ParticlesData_get_state(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+200);
  offset+=16+i0*8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set_state(ParticlesData obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+200);
  offset+=16+i0*8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp1_state(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+200);
  offset+=16+i0*8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  ArrNInt64 ParticlesData_getp_parent_particle_id(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+208);
  return (ArrNInt64)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len_parent_particle_id(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+208);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  int64_t ParticlesData_get_parent_particle_id(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+208);
  offset+=16+i0*8;
  return *( __global int64_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set_parent_particle_id(ParticlesData obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+208);
  offset+=16+i0*8;
  *( __global int64_t*)(( __global char*) obj+offset)=value;
}
   __global int64_t* ParticlesData_getp1_parent_particle_id(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+208);
  offset+=16+i0*8;
  return ( __global int64_t*)(( __global char*) obj+offset);
}
  ArrNUint32 ParticlesData_getp__rng_s1(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+216);
  return (ArrNUint32)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len__rng_s1(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+216);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  uint32_t ParticlesData_get__rng_s1(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+216);
  offset+=16+i0*4;
  return *( __global uint32_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__rng_s1(ParticlesData obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+216);
  offset+=16+i0*4;
  *( __global uint32_t*)(( __global char*) obj+offset)=value;
}
   __global uint32_t* ParticlesData_getp1__rng_s1(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+216);
  offset+=16+i0*4;
  return ( __global uint32_t*)(( __global char*) obj+offset);
}
  ArrNUint32 ParticlesData_getp__rng_s2(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+224);
  return (ArrNUint32)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len__rng_s2(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+224);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  uint32_t ParticlesData_get__rng_s2(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+224);
  offset+=16+i0*4;
  return *( __global uint32_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__rng_s2(ParticlesData obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+224);
  offset+=16+i0*4;
  *( __global uint32_t*)(( __global char*) obj+offset)=value;
}
   __global uint32_t* ParticlesData_getp1__rng_s2(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+224);
  offset+=16+i0*4;
  return ( __global uint32_t*)(( __global char*) obj+offset);
}
  ArrNUint32 ParticlesData_getp__rng_s3(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+232);
  return (ArrNUint32)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len__rng_s3(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+232);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  uint32_t ParticlesData_get__rng_s3(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+232);
  offset+=16+i0*4;
  return *( __global uint32_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__rng_s3(ParticlesData obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+232);
  offset+=16+i0*4;
  *( __global uint32_t*)(( __global char*) obj+offset)=value;
}
   __global uint32_t* ParticlesData_getp1__rng_s3(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+232);
  offset+=16+i0*4;
  return ( __global uint32_t*)(( __global char*) obj+offset);
}
  ArrNUint32 ParticlesData_getp__rng_s4(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+240);
  return (ArrNUint32)(( __global char*) obj+offset);
}
  int64_t ParticlesData_len__rng_s4(ParticlesData obj){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+240);
   __global int64_t* arr= ( __global int64_t*)(( __global char*) obj+offset);
  return arr[1];
}
  uint32_t ParticlesData_get__rng_s4(const ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+240);
  offset+=16+i0*4;
  return *( __global uint32_t*)(( __global char*) obj+offset);
}
  void ParticlesData_set__rng_s4(ParticlesData obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+240);
  offset+=16+i0*4;
  *( __global uint32_t*)(( __global char*) obj+offset)=value;
}
   __global uint32_t* ParticlesData_getp1__rng_s4(ParticlesData obj, int64_t i0){
  int64_t offset=0;
  offset+=*( __global int64_t*)(( __global char*) obj+offset+240);
  offset+=16+i0*4;
  return ( __global uint32_t*)(( __global char*) obj+offset);
}
#endif
// copyright ################################# //
// This file is part of the Xfields Package.   //
// Copyright (c) CERN, 2021.                   //
// ########################################### //

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif

#if !defined( TWO_OVER_SQRT_PI )
    #define TWO_OVER_SQRT_PI (1.128379167095512573896158903121545171688101258657997713688171443418)
#endif

#if !defined( SQRT_TWO )
    #define SQRT_TWO (1.414213562373095048801688724209698078569671875376948073176679738)
#endif

#if !defined( REAL_EPSILON )
    #define REAL_EPSILON 2.22044604925031e-16
#endif /* !defined( REAL_EPSILON ) */

// copyright ################################# //
// This file is part of the Xfields Package.   //
// Copyright (c) CERN, 2021.                   //
// ########################################### //

#ifndef CENTRAL_DIFF_H
#define CENTRAL_DIFF_H

 __kernel 
void central_diff(
	      const int     nelem,
	      const int     row_size,
	      const int     stride_in_dbl,
	      const double  factor,
 __global   const int8_t* matrix_buffer,
              const int64_t matrix_offset,
 __global         int8_t* res_buffer,
                    int64_t res_offset
              ){

    __global  const double* matrix = 
	           ( __global  double*) (matrix_buffer + matrix_offset); 
    __global        double*  res = 
	           ( __global  double*) (res_buffer + res_offset); 

int ii; //autovectorized

ii=get_global_id(0); //autovectorized

      if (ii-stride_in_dbl>=0 && ii+stride_in_dbl<nelem){
         res[ii] = factor * (matrix[ii+stride_in_dbl]
			   - matrix[ii-stride_in_dbl]);
      } 
      int place_in_row = (ii/stride_in_dbl)%row_size;
      if (place_in_row==0 || place_in_row==row_size-1){
         res[ii] = 0;
      }
//end autovectorized


}

#endif

// copyright ################################# //
// This file is part of the Xfields Package.   //
// Copyright (c) CERN, 2021.                   //
// ########################################### //

#ifndef XFIELDS_LINEAR_INTERPOLATORS_H
#define XFIELDS_LINEAR_INTERPOLATORS_H

typedef struct{
    int64_t ix;
    int64_t iy;
    int64_t iz;
    int64_t nx;
    int64_t ny;
    int64_t nz;
    double w000;
    double w100;
    double w010;
    double w110;
    double w001;
    double w101;
    double w011;
    double w111;
}IndicesAndWeights;


 
IndicesAndWeights TriLinearInterpolatedFieldMap_compute_indeces_and_weights(
	TriLinearInterpolatedFieldMapData fmap,
	double x, double y, double z){

	IndicesAndWeights iw;

	const double dx = TriLinearInterpolatedFieldMapData_get_dx(fmap);
	const double dy = TriLinearInterpolatedFieldMapData_get_dy(fmap);
	const double dz = TriLinearInterpolatedFieldMapData_get_dz(fmap);
	const double x0 = TriLinearInterpolatedFieldMapData_get_x_min(fmap);
	const double y0 = TriLinearInterpolatedFieldMapData_get_y_min(fmap);
	const double z0 = TriLinearInterpolatedFieldMapData_get_z_min(fmap);
	const int64_t nx = TriLinearInterpolatedFieldMapData_get_nx(fmap);
	const int64_t ny = TriLinearInterpolatedFieldMapData_get_ny(fmap);
	const int64_t nz = TriLinearInterpolatedFieldMapData_get_nz(fmap);

    	iw.nx = nx;
    	iw.ny = ny;
    	iw.nz = nz;

    	// indices
    	iw.ix = floor((x - x0) / dx);
    	iw.iy = floor((y - y0) / dy);
    	iw.iz = floor((z - z0) / dz);

	
    	if (iw.ix >= 0 && iw.ix < nx - 1 && iw.iy >= 0 && iw.iy < ny - 1
	    	    && iw.iz >= 0 && iw.iz < nz - 1){

    	    // distances
    	    const double dxi = x - (x0 + iw.ix * dx);
    	    const double dyi = y - (y0 + iw.iy * dy);
    	    const double dzi = z - (z0 + iw.iz * dz);
	    
    	    // weights
    	    iw.w000 = (1.-dxi/dx) * (1.-dyi/dy) * (1.-dzi/dz);
    	    iw.w100 = (dxi/dx)    * (1.-dyi/dy) * (1.-dzi/dz);
    	    iw.w010 = (1.-dxi/dx) * (dyi/dy)    * (1.-dzi/dz);
    	    iw.w110 = (dxi/dx)    * (dyi/dy)    * (1.-dzi/dz);
    	    iw.w001 = (1.-dxi/dx) * (1.-dyi/dy) * (dzi/dz);
    	    iw.w101 = (dxi/dx)    * (1.-dyi/dy) * (dzi/dz);
    	    iw.w011 = (1.-dxi/dx) * (dyi/dy)    * (dzi/dz);
    	    iw.w111 = (dxi/dx)    * (dyi/dy)    * (dzi/dz);
	}
	else{
            iw.ix = -999; 
            iw.iy = -999; 
            iw.iz = -999; 
	}
	return iw;

}	

 
double TriLinearInterpolatedFieldMap_interpolate_3d_map_scalar(
	 __global  const double* map,
	   const IndicesAndWeights iw){
	
    double val;

    if (iw.ix < 0){
	 val = 0.;
    }
    else{
	val = 
    	       iw.w000 * map[iw.ix   + (iw.iy  ) * iw.nx + (iw.iz  ) * iw.nx * iw.ny]
    	     + iw.w100 * map[iw.ix+1 + (iw.iy  ) * iw.nx + (iw.iz  ) * iw.nx * iw.ny]
    	     + iw.w010 * map[iw.ix+  + (iw.iy+1) * iw.nx + (iw.iz  ) * iw.nx * iw.ny]
    	     + iw.w110 * map[iw.ix+1 + (iw.iy+1) * iw.nx + (iw.iz  ) * iw.nx * iw.ny]
    	     + iw.w001 * map[iw.ix   + (iw.iy  ) * iw.nx + (iw.iz+1) * iw.nx * iw.ny]
    	     + iw.w101 * map[iw.ix+1 + (iw.iy  ) * iw.nx + (iw.iz+1) * iw.nx * iw.ny]
    	     + iw.w011 * map[iw.ix+  + (iw.iy+1) * iw.nx + (iw.iz+1) * iw.nx * iw.ny]
    	     + iw.w111 * map[iw.ix+1 + (iw.iy+1) * iw.nx + (iw.iz+1) * iw.nx * iw.ny];
    }

    return val;
}

 __kernel 
void TriLinearInterpolatedFieldMap_interpolate_3d_map_vector(
    TriLinearInterpolatedFieldMapData  fmap,
                        const int64_t  n_points,
            __global  const double*  x,
            __global  const double*  y,
            __global  const double*  z,
                        const int64_t  n_quantities,
            __global  const int8_t*  buffer_mesh_quantities,
            __global  const int64_t* offsets_mesh_quantities,
            __global        double*  particles_quantities) {

//    #pragma omp parallel for //only_for_context cpu_openmp 
int pidx; //autovectorized

pidx=get_global_id(0); //autovectorized


	const IndicesAndWeights iw = 
		TriLinearInterpolatedFieldMap_compute_indeces_and_weights(
	                                      fmap, x[pidx], y[pidx], z[pidx]);
    	for (int iq=0; iq<n_quantities; iq++){
	    particles_quantities[iq*n_points + pidx] = 
		TriLinearInterpolatedFieldMap_interpolate_3d_map_scalar(
	           ( __global  double*)(buffer_mesh_quantities + offsets_mesh_quantities[iq]),
		   iw);
	}
//end autovectorized

}
#endif

// copyright ################################# //
// This file is part of the Xfields Package.   //
// Copyright (c) CERN, 2021.                   //
// ########################################### //

#ifndef XFIELDS_CHARGE_DEPOSITION_H
#define XFIELDS_CHARGE_DEPOSITION_H


//from file: atomicadd.clh

#ifndef _ATOMICADD_CLH_
#define _ATOMICADD_CLH_


#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
inline void atomicAdd(volatile __global double *addr, double val)
{
	union {
		long u64;
		double f64;
	} next, expected, current;
	current.f64 = *addr;
	do {
		expected.f64 = current.f64;
		next.f64 = expected.f64 + val;
		current.u64 = atom_cmpxchg(
			(volatile __global long *)addr,
		        (long) expected.u64,
			(long) next.u64);
	} while( current.u64 != expected.u64 );
}
#endif

//end file: atomicadd.clh





  void p2m_rectmesh3d_one_particle(
        // INPUTS:
        const double x, 
	const double y, 
	const double z,
	  // particle weight
	const double pwei,
          // mesh origin
        const double x0, const double y0, const double z0,
          // mesh distances per cell
        const double dx, const double dy, const double dz,
          // mesh dimension (number of cells)
        const int nx, const int ny, const int nz,
        // OUTPUTS:
         __global  double *grid1d
) {

    double vol_m1 = 1/(dx*dy*dz);

    // indices
    int jx = floor((x - x0) / dx);
    int ix = floor((y - y0) / dy);
    int kx = floor((z - z0) / dz);

    // distances
    double dxi = x - (x0 + jx * dx);
    double dyi = y - (y0 + ix * dy);
    double dzi = z - (z0 + kx * dz);

    // weights
    double wijk =    pwei * vol_m1 * (1.-dxi/dx) * (1.-dyi/dy) * (1.-dzi/dz);
    double wi1jk =   pwei * vol_m1 * (1.-dxi/dx) * (dyi/dy)    * (1.-dzi/dz);
    double wij1k =   pwei * vol_m1 * (dxi/dx)    * (1.-dyi/dy) * (1.-dzi/dz);
    double wi1j1k =  pwei * vol_m1 * (dxi/dx)    * (dyi/dy)    * (1.-dzi/dz);
    double wijk1 =   pwei * vol_m1 * (1.-dxi/dx) * (1.-dyi/dy) * (dzi/dz);
    double wi1jk1 =  pwei * vol_m1 * (1.-dxi/dx) * (dyi/dy)    * (dzi/dz);
    double wij1k1 =  pwei * vol_m1 * (dxi/dx)    * (1.-dyi/dy) * (dzi/dz);
    double wi1j1k1 = pwei * vol_m1 * (dxi/dx)    * (dyi/dy)    * (dzi/dz);

    if (jx >= 0 && jx < nx - 1 && ix >= 0 && ix < ny - 1
        	    && kx >= 0 && kx < nz - 1)
    {
        atomicAdd(&grid1d[jx   + ix*nx     + kx*nx*ny],     wijk);
        atomicAdd(&grid1d[jx+1 + ix*nx     + kx*nx*ny],     wij1k);
        atomicAdd(&grid1d[jx   + (ix+1)*nx + kx*nx*ny],     wi1jk);
        atomicAdd(&grid1d[jx+1 + (ix+1)*nx + kx*nx*ny],     wi1j1k);
        atomicAdd(&grid1d[jx   + ix*nx     + (kx+1)*nx*ny], wijk1);
        atomicAdd(&grid1d[jx+1 + ix*nx     + (kx+1)*nx*ny], wij1k1);
        atomicAdd(&grid1d[jx   + (ix+1)*nx + (kx+1)*nx*ny], wi1jk1);
        atomicAdd(&grid1d[jx+1 + (ix+1)*nx + (kx+1)*nx*ny], wi1j1k1);
    }

}


 __kernel  void p2m_rectmesh3d(
        // INPUTS:
          // length of x, y, z arrays
        const int nparticles,
          // particle positions
         __global  const double* x, 
	 __global  const double* y, 
	 __global  const double* z,
	  // particle weights and stat flags
	 __global  const double* part_weights,
	 __global  const int64_t* part_state,
          // mesh origin
        const double x0, const double y0, const double z0,
          // mesh distances per cell
        const double dx, const double dy, const double dz,
          // mesh dimension (number of cells)
        const int nx, const int ny, const int nz,
        // OUTPUTS:
         __global  int8_t*  grid1d_buffer,
	             int64_t  grid1d_offset){

     __global  double* grid1d = 
		( __global  double*)(grid1d_buffer + grid1d_offset);

//    #pragma omp parallel for //only_for_context cpu_openmp 
int pidx; //autovectorized

pidx=get_global_id(0); //autovectorized

        if (part_state[pidx] > 0){
    	    double pwei = part_weights[pidx];

            p2m_rectmesh3d_one_particle(x[pidx], y[pidx], z[pidx], pwei,
                                        x0, y0, z0, dx, dy, dz, nx, ny, nz,
                                        grid1d);
	}
//end autovectorized

}

 __kernel  void p2m_rectmesh3d_xparticles(
        // INPUTS:
          // length of x, y, z arrays
        const int nparticles,
	ParticlesData particles,
          // mesh origin
        const double x0, const double y0, const double z0,
          // mesh distances per cell
        const double dx, const double dy, const double dz,
          // mesh dimension (number of cells)
        const int nx, const int ny, const int nz,
        // OUTPUTS:
         __global  int8_t*  grid1d_buffer,
	             int64_t  grid1d_offset){

     __global  double* grid1d = 
    	( __global  double*)(grid1d_buffer + grid1d_offset);
    
     __global  const double* x = ParticlesData_getp1_x(particles, 0); 
     __global  const double* y = ParticlesData_getp1_y(particles, 0); 
     __global  const double* z = ParticlesData_getp1_zeta(particles, 0);
     __global  const double* part_weights = ParticlesData_getp1_weight(
    		                                             particles, 0);
     __global  const int64_t* part_state = ParticlesData_getp1_state(
    		                                             particles, 0);
    // TODO I am forgetting about charge_ratio and mass_ratio
    const double q0_coulomb = QELEM * ParticlesData_get_q0(particles);

//    #pragma omp parallel for //only_for_context cpu_openmp 
int pidx; //autovectorized

pidx=get_global_id(0); //autovectorized

        if (part_state[pidx] > 0){
    	    double pwei = part_weights[pidx] * q0_coulomb;

            p2m_rectmesh3d_one_particle(x[pidx], y[pidx], z[pidx], pwei,
                                        x0, y0, z0, dx, dy, dz, nx, ny, nz,
                                        grid1d);
	}
//end autovectorized


}
#endif