#pragma once

#include <float.h>
#include <math.h>
#include <vector>

// NOTE:
// Ported from AntTweakBar
struct ImOrient
{
    typedef unsigned int color32;

    const ImU32 COLOR32_BLACK = 0xff000000;   // Black 
    const ImU32 COLOR32_WHITE = 0xffffffff;   // White 
    const ImU32 COLOR32_ZERO = 0x00000000;    // Zero 
    const ImU32 COLOR32_RED = 0xffff0000;     // Red 
    const ImU32 COLOR32_GREEN = 0xff00ff00;   // Green 
    const ImU32 COLOR32_BLUE = 0xff0000ff;    // Blue 

    enum Cull { CULL_NONE, CULL_CW, CULL_CCW };

    IMGUI_API ImOrient();
    IMGUI_API bool Orient(char* label);

    // TODO: Clean these up; not all necessary
    double               Qx, Qy, Qz, Qs;    // Quat value
    double               Vx, Vy, Vz, Angle; // Not used
    double               Dx, Dy, Dz;        // Dir value set when used as a direction
    bool                 m_AAMode;          // Axis & angle mode -> disabled
    bool                 m_ShowVal;         // Display values
    bool                 m_IsFloat;         // Quat/Dir uses floats
    bool                 m_IsDir;           // Mapped to a dir vector instead of a quat
    double               m_Dir[3];          // If not zero, display one direction vector
    ImU32               m_DirColor;        // Direction vector color
    float                m_Permute[3][3];   // Permute frame axis
    bool                 m_Highlighted;
    bool                 m_Rotating;
    double               m_OrigQuat[4];
    float                m_OrigX, m_OrigY;
    double               m_PrevX, m_PrevY;

    // For the geometry
    enum EArrowParts { ARROW_CONE, ARROW_CONE_CAP, ARROW_CYL, ARROW_CYL_CAP };
    static std::vector<float>   s_SphTri;
    static std::vector<ImU32> s_SphCol;
    static std::vector<int>     s_SphTriProj;
    static std::vector<ImU32> s_SphColLight;
    static std::vector<float>   s_ArrowTri[4];
    static std::vector<int>     s_ArrowTriProj[4];
    static std::vector<float>   s_ArrowNorm[4];
    static std::vector<ImU32> s_ArrowColLight[4];
    static void CreateSphere();
    static void CreateArrow();

    IMGUI_API void ConvertToAxisAngle();
    IMGUI_API void ConvertFromAxisAngle();
    IMGUI_API void ApplyQuat(float *outX, float *outY, float *outZ, float x, float y, float z, float qx, float qy, float qz, float qs);
    IMGUI_API void QuatFromDir(double *outQx, double *outQy, double *outQz, double *outQs, double dx, double dy, double dz);
    IMGUI_API void Permute(float *outX, float *outY, float *outZ, float x, float y, float z);
    IMGUI_API void PermuteInv(float *outX, float *outY, float *outZ, float x, float y, float z);
    IMGUI_API void Permute(double *outX, double *outY, double *outZ, double x, double y, double z);
    IMGUI_API void PermuteInv(double *outX, double *outY, double *outZ, double x, double y, double z);

    inline double DegToRad(double degree) { return degree * (M_PI / 180.0); }
    inline double RadToDeg(double radian) { return radian * (180.0 / M_PI); }

    IMGUI_API static ImU32 ColorBlend(ImU32 _Color1, ImU32 _Color2, float _S);
    IMGUI_API static void QuatMult(double *out, const double *q1, const double *q2);
    IMGUI_API static void QuatFromAxisAngle(double *out, const double *axis, double angle);
    IMGUI_API static void Vec3Cross(double *out, const double *a, const double *b);
    IMGUI_API static double Vec3Dot(const double *a, const double *b);
    IMGUI_API static void Vec3RotY(float *x, float *y, float *z);
    IMGUI_API static void Vec3RotZ(float *x, float *y, float *z);
    IMGUI_API static float QuatD(int w, int h);
    IMGUI_API static int QuatPX(float x, int w, int h);
    IMGUI_API static int QuatPY(float y, int w, int h);
    IMGUI_API static float QuatIX(int x, int w, int h);
    IMGUI_API static float QuatIY(int y, int w, int h);
     
    /*
    static void IMGUI_API CopyVarFromExtCB(void *_VarValue, const void *_ExtValue, unsigned int _ExtMemberIndex, void *_ClientData);
    static void IMGUI_API CopyVarToExtCB(const void *_VarValue, void *_ExtValue, unsigned int _ExtMemberIndex, void *_ClientData);
    static void IMGUI_API SummaryCB(char *_SummaryString, size_t _SummaryMaxLength, const void *_ExtValue, void *_ClientData);
    static void IMGUI_API DrawCB(int _W, int _H, void *_ExtValue, void *_ClientData);
    static bool IMGUI_API MouseMotionCB(int _MouseX, int _MouseY, int _W, int _H, void *_StructExtValue, void *_ClientData);
    static bool IMGUI_API MouseButtonCB(TwMouseButtonID _Button, bool _Pressed, int _MouseX, int _MouseY, int _W, int _H, void *_StructExtValue, void *_ClientData);
    static void IMGUI_API MouseLeaveCB(void *_StructExtValue, void *_ClientData);
    */
};

template <typename _T> inline const _T& TClamp(const _T& _X, const _T& _Limit1, const _T& _Limit2)
{
    if (_Limit1 < _Limit2)
        return (_X <= _Limit1) ? _Limit1 : ((_X >= _Limit2) ? _Limit2 : _X);
    else
        return (_X <= _Limit2) ? _Limit2 : ((_X >= _Limit1) ? _Limit1 : _X);
}

