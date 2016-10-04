#include "imgui.h"
#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui_internal.h" // ImSaturate
#include "imgui_orient.h"

ImVector<float> ImOrient::s_SphTri;
ImVector<ImU32> ImOrient::s_SphCol;
ImVector<int>   ImOrient::s_SphTriProj;
ImVector<ImU32> ImOrient::s_SphColLight;
ImVector<float> ImOrient::s_ArrowTri[4];
ImVector<int>   ImOrient::s_ArrowTriProj[4];
ImVector<float> ImOrient::s_ArrowNorm[4];
ImVector<ImU32> ImOrient::s_ArrowColLight[4];

bool ImOrient::Orient(char* label)
{
    //    ImGuiIO& io = ImGui::GetIO();
    ImGuiStyle& style = ImGui::GetStyle();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();

    if (ImOrient::s_SphTri.empty())
    {
        ImOrient::CreateArrow();
        ImOrient::CreateSphere();
    }

    ImGui::PushID(label);
    ImGui::BeginGroup();

    bool value_changed = false;

    if (m_AAMode)
    {
        ImGui::Text("V={%.2f,%.2f,%.2f} A=%.0f%c", Vx, Vy, Vz, Angle, 176);
    }
    else if (m_IsDir)
    {
        ImGui::Text("V={%.2f,%.2f,%.2f}", m_Dir[0], m_Dir[1], m_Dir[2]);
    }
    else
    {
        ImGui::Text("Q={x:%.2f,y:%.2f,z:%.2f,s:%.2f}", Qx, Qy, Qz, Qs);
    }

    ImVec2 orient_pos = ImGui::GetCursorScreenPos();
    float sv_orient_size = ImGui::CalcItemWidth() / 3.0f;
    ImGui::InvisibleButton("orient", ImVec2(sv_orient_size, sv_orient_size));

    draw_list->AddRectFilled(orient_pos, orient_pos + ImVec2(sv_orient_size, sv_orient_size), ImColor(style.Colors[ImGuiCol_FrameBg]), style.FrameRounding);

    float w = sv_orient_size;
    float h = sv_orient_size;

    double normDir = sqrt(m_Dir[0] * m_Dir[0] + m_Dir[1] * m_Dir[1] + m_Dir[2] * m_Dir[2]);
    bool drawDir = m_IsDir || (normDir > DBL_EPSILON);

    ImVec2 inner_pos = orient_pos;
    float inner_size = w;
    if (drawDir)
    {
        inner_size = sv_orient_size;
    }
    else
    {
        inner_pos.x += sv_orient_size * .25f * .5f;
        inner_pos.y += sv_orient_size * .25f * .5f;
        inner_size *= .75f;
    }
    float x, y, z, nx, ny, nz, kx, ky, kz, qx, qy, qz, qs;
    int i, j, k, l, m;

    // normalize quaternion
    float qn = (float)sqrt(Qs*Qs + Qx*Qx + Qy*Qy + Qz*Qz);
    if (qn > FLT_EPSILON)
    {
        qx = (float)Qx / qn;
        qy = (float)Qy / qn;
        qz = (float)Qz / qn;
        qs = (float)Qs / qn;
    }
    else
    {
        qx = qy = qz = 0;
        qs = 1;
    }

    ImColor alpha(1.0f, 1.0f, 1.0f, m_Highlighted ? 1.0f : 0.74f);

    // check if frame is right-handed
    Permute(&kx, &ky, &kz, 1, 0, 0);
    double px[3] = { (double)kx, (double)ky, (double)kz };
    Permute(&kx, &ky, &kz, 0, 1, 0);
    double py[3] = { (double)kx, (double)ky, (double)kz };
    Permute(&kx, &ky, &kz, 0, 0, 1);
    double pz[3] = { (double)kx, (double)ky, (double)kz };
    double ez[3];
    Vec3Cross(ez, px, py);

    bool AA = style.AntiAliasedShapes;
    style.AntiAliasedShapes = false;
    bool frameRightHanded = (ez[0] * pz[0] + ez[1] * pz[1] + ez[2] * pz[2] >= 0);
    float cullDir = frameRightHanded ? 1.0f : -1.0f;

    if (drawDir)
    {
        float dir[] = { (float)m_Dir[0], (float)m_Dir[1], (float)m_Dir[2] };
        if (normDir < DBL_EPSILON)
        {
            normDir = 1;
            dir[0] = 1;
        }
        kx = dir[0]; ky = dir[1]; kz = dir[2];
        double rotDirAxis[3] = { 0, -kz, ky };
        if (rotDirAxis[0] * rotDirAxis[0] + rotDirAxis[1] * rotDirAxis[1] + rotDirAxis[2] * rotDirAxis[2] < (DBL_EPSILON*DBL_EPSILON))
        {
            rotDirAxis[0] = rotDirAxis[1] = 0;
            rotDirAxis[2] = 1;
        }
        double rotDirAngle = acos(kx / normDir);
        double rotDirQuat[4];
        QuatFromAxisAngle(rotDirQuat, rotDirAxis, rotDirAngle);

        kx = 1; ky = 0; kz = 0;
        ApplyQuat(&kx, &ky, &kz, kx, ky, kz, (float)rotDirQuat[0], (float)rotDirQuat[1], (float)rotDirQuat[2], (float)rotDirQuat[3]);
        ApplyQuat(&kx, &ky, &kz, kx, ky, kz, qx, qy, qz, qs);
        for (k = 0; k < 4; ++k) // 4 parts of the arrow
        {
            // draw order
            Permute(&x, &y, &z, kx, ky, kz);
            j = (z > 0) ? 3 - k : k;

            assert(s_ArrowTriProj[j].size() == 2 * (s_ArrowTri[j].size() / 3) && s_ArrowColLight[j].size() == s_ArrowTri[j].size() / 3 && s_ArrowNorm[j].size() == s_ArrowTri[j].size());
            const int ntri = (int)s_ArrowTri[j].size() / 3;
            const float *tri = &(s_ArrowTri[j][0]);
            const float *norm = &(s_ArrowNorm[j][0]);
            int *triProj = &(s_ArrowTriProj[j][0]);
            ImU32 *colLight = &(s_ArrowColLight[j][0]);
            for (i = 0; i < ntri; ++i)
            {
                x = tri[3 * i + 0]; y = tri[3 * i + 1]; z = tri[3 * i + 2];
                nx = norm[3 * i + 0]; ny = norm[3 * i + 1]; nz = norm[3 * i + 2];
                if (x > 0)
                    x = 2.5f*x - 2.0f;
                else
                    x += 0.2f;
                y *= 1.5f;
                z *= 1.5f;
                ApplyQuat(&x, &y, &z, x, y, z, (float)rotDirQuat[0], (float)rotDirQuat[1], (float)rotDirQuat[2], (float)rotDirQuat[3]);
                ApplyQuat(&x, &y, &z, x, y, z, qx, qy, qz, qs);
                Permute(&x, &y, &z, x, y, z);
                ApplyQuat(&nx, &ny, &nz, nx, ny, nz, (float)rotDirQuat[0], (float)rotDirQuat[1], (float)rotDirQuat[2], (float)rotDirQuat[3]);
                ApplyQuat(&nx, &ny, &nz, nx, ny, nz, qx, qy, qz, qs);
                Permute(&nx, &ny, &nz, nx, ny, nz);
                triProj[2 * i + 0] = QuatPX(x, int(w), int(h));
                triProj[2 * i + 1] = QuatPY(y, int(w), int(h));
                ImU32 col = (m_DirColor | 0xff000000) & alpha;
                colLight[i] = ColorBlend(0xff000000, col, fabsf(TClamp(nz, -1.0f, 1.0f)));
            }

            DrawTriangles(draw_list, inner_pos, triProj, colLight, ntri, cullDir);
        }
    }
    else
    {
        // draw arrows & sphere
        const float SPH_RADIUS = 0.75f;
        for (m = 0; m < 2; ++m)  // m=0: back, m=1: front
        {
            for (l = 0; l < 3; ++l)  // draw 3 arrows
            {
                kx = 1; ky = 0; kz = 0;
                if (l == 1)
                    Vec3RotZ(&kx, &ky, &kz);
                else if (l == 2)
                    Vec3RotY(&kx, &ky, &kz);
                ImOrient::ApplyQuat(&kx, &ky, &kz, kx, ky, kz, qx, qy, qz, qs);
                for (k = 0; k < 4; ++k) // 4 parts of the arrow
                {
                    // draw order
                    Permute(&x, &y, &z, kx, ky, kz);
                    j = (z > 0) ? 3 - k : k;

                    bool cone = true;
                    if ((m == 0 && z > 0) || (m == 1 && z <= 0))
                    {
                        if (j == ImOrient::ARROW_CONE || j == ImOrient::ARROW_CONE_CAP) // do not draw cone
                            continue;
                        else
                            cone = false;
                    }
                    assert(ImOrient::s_ArrowTriProj[j].size() == 2 * (ImOrient::s_ArrowTri[j].size() / 3) && ImOrient::s_ArrowColLight[j].size() == ImOrient::s_ArrowTri[j].size() / 3 && ImOrient::s_ArrowNorm[j].size() == ImOrient::s_ArrowTri[j].size());
                    const int ntri = (int)ImOrient::s_ArrowTri[j].size() / 3;
                    const float *tri = &(ImOrient::s_ArrowTri[j][0]);
                    const float *norm = &(ImOrient::s_ArrowNorm[j][0]);
                    int *triProj = &(ImOrient::s_ArrowTriProj[j][0]);
                    ImU32 *colLight = &(ImOrient::s_ArrowColLight[j][0]);
                    for (i = 0; i < ntri; ++i)
                    {
                        x = tri[3 * i + 0]; y = tri[3 * i + 1]; z = tri[3 * i + 2];
                        if (cone && x <= 0)
                            x = SPH_RADIUS;
                        else if (!cone && x > 0)
                            x = -SPH_RADIUS;
                        nx = norm[3 * i + 0]; ny = norm[3 * i + 1]; nz = norm[3 * i + 2];
                        if (l == 1)
                        {
                            Vec3RotZ(&x, &y, &z);
                            Vec3RotZ(&nx, &ny, &nz);
                        }
                        else if (l == 2)
                        {
                            Vec3RotY(&x, &y, &z);
                            Vec3RotY(&nx, &ny, &nz);
                        }
                        ImOrient::ApplyQuat(&x, &y, &z, x, y, z, qx, qy, qz, qs);
                        Permute(&x, &y, &z, x, y, z);
                        ImOrient::ApplyQuat(&nx, &ny, &nz, nx, ny, nz, qx, qy, qz, qs);
                        Permute(&nx, &ny, &nz, nx, ny, nz);
                        triProj[2 * i + 0] = QuatPX(x, int(inner_size), int(inner_size));
                        triProj[2 * i + 1] = QuatPY(y, int(inner_size), int(inner_size));
                        float fade = (m == 0 && z < 0) ? TClamp(2.0f*z*z, 0.0f, 1.0f) : 0;
                        float alphaFade = 1.0f;
                        alphaFade = alpha.Value.w;
                        alphaFade *= (1.0f - fade);
                        ImColor alphaFadeCol(1.0f, 1.0f, 1.0f, alphaFade);
                        ImU32 col = (l == 0) ? 0xffff0000 : ((l == 1) ? 0xff00ff00 : 0xff0000ff);
                        colLight[i] = ColorBlend(0xff000000, col, fabsf(TClamp(nz, -1.0f, 1.0f))) & ImU32(alphaFadeCol);
                    }

                    DrawTriangles(draw_list, inner_pos, triProj, colLight, ntri, cullDir);
                }
            }

            if (m == 0)
            {
                const float *tri = &(ImOrient::s_SphTri[0]);
                int *triProj = &(ImOrient::s_SphTriProj[0]);
                const ImU32 *col = &(ImOrient::s_SphCol[0]);
                ImU32 *colLight = &(ImOrient::s_SphColLight[0]);
                const int ntri = (int)ImOrient::s_SphTri.size() / 3;
                for (i = 0; i < ntri; ++i)   // draw sphere
                {
                    x = SPH_RADIUS*tri[3 * i + 0]; y = SPH_RADIUS*tri[3 * i + 1]; z = SPH_RADIUS*tri[3 * i + 2];
                    ImOrient::ApplyQuat(&x, &y, &z, x, y, z, qx, qy, qz, qs);
                    Permute(&x, &y, &z, x, y, z);
                    triProj[2 * i + 0] = QuatPX(x, int(inner_size), int(inner_size));
                    triProj[2 * i + 1] = QuatPY(y, int(inner_size), int(inner_size));
                    colLight[i] = ColorBlend(0xff000000, col[i], fabsf(TClamp(z / SPH_RADIUS, -1.0f, 1.0f)))& ImU32(alpha);
                }

                DrawTriangles(draw_list, inner_pos, triProj, colLight, ntri, cullDir);
            }
        }

        // draw x
        draw_list->AddLine(orient_pos + ImVec2(w - 12, h - 36), orient_pos + ImVec2(w - 12 + 5, h - 36 + 5), 0xffc00000);
        draw_list->AddLine(orient_pos + ImVec2(w - 12 + 5, h - 36), orient_pos + ImVec2(w - 12, h - 36 + 5), 0xffc00000);
        // draw y
        draw_list->AddLine(orient_pos + ImVec2(w - 12, h - 25), orient_pos + ImVec2(w - 12 + 3, h - 25 + 4), 0xff00c000);
        draw_list->AddLine(orient_pos + ImVec2(w - 12 + 5, h - 25), orient_pos + ImVec2(w - 12, h - 25 + 7), 0xff00c000);
        // draw z
        draw_list->AddLine(orient_pos + ImVec2(w - 12, h - 12), orient_pos + ImVec2(w - 12 + 5, h - 12), 0xff0000c0);
        draw_list->AddLine(orient_pos + ImVec2(w - 12, h - 12 + 5), orient_pos + ImVec2(w - 12 + 5, h - 12 + 5), 0xff0000c0);
        draw_list->AddLine(orient_pos + ImVec2(w - 12, h - 12 + 5), orient_pos + ImVec2(w - 12 + 5, h - 12), 0xff0000c0);
    }

    ImGui::EndGroup();
    ImGui::PopID();


    style.AntiAliasedShapes = AA;
    return value_changed;
}

void ImOrient::DrawTriangles(ImDrawList* draw_list, ImVec2 offset, int* triProj, ImU32* colLight, int numVertices, float cullDir)
{
    const ImVec2 uv = GImGui->FontTexUvWhitePixel;
    draw_list->PrimReserve(numVertices, numVertices); // num vert/indices 
    for (int ii = 0; ii < numVertices / 3; ii++)
    {
        ImVec2 v1 = offset + ImVec2(float(triProj[ii * 6 + 0]), float(triProj[ii * 6 + 1]));
        ImVec2 v2 = offset + ImVec2(float(triProj[ii * 6 + 2]), float(triProj[ii * 6 + 3]));
        ImVec2 v3 = offset + ImVec2(float(triProj[ii * 6 + 4]), float(triProj[ii * 6 + 5]));

        // 2D cross product to do culling
        ImVec2 d1 = Vec2Subtract(v2, v1);
        ImVec2 d2 = Vec2Subtract(v3, v1);
        float c = Vec2Cross(d1, d2) * cullDir;
        if (c > 0.0f)
        {
            v2 = v1;
            v3 = v1;
        }

        draw_list->PrimWriteIdx(ImDrawIdx(draw_list->_VtxCurrentIdx));
        draw_list->PrimWriteIdx(ImDrawIdx(draw_list->_VtxCurrentIdx + 1));
        draw_list->PrimWriteIdx(ImDrawIdx(draw_list->_VtxCurrentIdx + 2));
        draw_list->PrimWriteVtx(v1, uv, colLight[ii * 3]);
        draw_list->PrimWriteVtx(v2, uv, colLight[ii * 3 + 1]);
        draw_list->PrimWriteVtx(v3, uv, colLight[ii * 3 + 2]);
    }
}
            
void ImOrient::CreateSphere()
{
    const int SUBDIV = 7;
    s_SphTri.clear();
    s_SphCol.clear();

    const float A[8 * 3] = { 1,0,0, 0,0,-1, -1,0,0, 0,0,1,   0,0,1,  1,0,0,  0,0,-1, -1,0,0 };
    const float B[8 * 3] = { 0,1,0, 0,1,0,  0,1,0,  0,1,0,   0,-1,0, 0,-1,0, 0,-1,0, 0,-1,0 };
    const float C[8 * 3] = { 0,0,1, 1,0,0,  0,0,-1, -1,0,0,  1,0,0,  0,0,-1, -1,0,0, 0,0,1 };
    const ImU32 COL_A[8] = { 0xffffffff, 0xff40ffFF, 0xff40ff40, 0xffffff40,  0xffff40ff, 0xff4040FF, 0xff404040, 0xffFF4040 };
    const ImU32 COL_B[8] = { 0xffffffff, 0xff40ffFF, 0xff40ff40, 0xffffff40,  0xffff40ff, 0xff4040FF, 0xff404040, 0xffFF4040 };
    const ImU32 COL_C[8] = { 0xffffffff, 0xff40ffFF, 0xff40ff40, 0xffffff40,  0xffff40ff, 0xff4040FF, 0xff404040, 0xffFF4040 };

    int i, j, k, l;
    float xa, ya, za, xb, yb, zb, xc, yc, zc, x, y, z, norm, u[3], v[3];
    ImU32 col;
    for (i = 0; i < 8; ++i)
    {
        xa = A[3 * i + 0]; ya = A[3 * i + 1]; za = A[3 * i + 2];
        xb = B[3 * i + 0]; yb = B[3 * i + 1]; zb = B[3 * i + 2];
        xc = C[3 * i + 0]; yc = C[3 * i + 1]; zc = C[3 * i + 2];
        for (j = 0; j <= SUBDIV; ++j)
            for (k = 0; k <= 2 * (SUBDIV - j); ++k)
            {
                if (k % 2 == 0)
                {
                    u[0] = ((float)j) / (SUBDIV + 1);
                    v[0] = ((float)(k / 2)) / (SUBDIV + 1);
                    u[1] = ((float)(j + 1)) / (SUBDIV + 1);
                    v[1] = ((float)(k / 2)) / (SUBDIV + 1);
                    u[2] = ((float)j) / (SUBDIV + 1);
                    v[2] = ((float)(k / 2 + 1)) / (SUBDIV + 1);
                }
                else
                {
                    u[0] = ((float)j) / (SUBDIV + 1);
                    v[0] = ((float)(k / 2 + 1)) / (SUBDIV + 1);
                    u[1] = ((float)(j + 1)) / (SUBDIV + 1);
                    v[1] = ((float)(k / 2)) / (SUBDIV + 1);
                    u[2] = ((float)(j + 1)) / (SUBDIV + 1);
                    v[2] = ((float)(k / 2 + 1)) / (SUBDIV + 1);
                }

                for (l = 0; l < 3; ++l)
                {
                    x = (1.0f - u[l] - v[l])*xa + u[l] * xb + v[l] * xc;
                    y = (1.0f - u[l] - v[l])*ya + u[l] * yb + v[l] * yc;
                    z = (1.0f - u[l] - v[l])*za + u[l] * zb + v[l] * zc;
                    norm = sqrtf(x*x + y*y + z*z);
                    x /= norm; y /= norm; z /= norm;
                    s_SphTri.push_back(x); s_SphTri.push_back(y); s_SphTri.push_back(z);
                    if (u[l] + v[l] > FLT_EPSILON)
                        col = ColorBlend(COL_A[i], ColorBlend(COL_B[i], COL_C[i], v[l] / (u[l] + v[l])), u[l] + v[l]);
                    else
                        col = COL_A[i];
                    //if( (j==0 && k==0) || (j==0 && k==2*SUBDIV) || (j==SUBDIV && k==0) )
                    //  col = 0xffff0000;
                    s_SphCol.push_back(col);
                }
            }
    }
    s_SphTriProj.clear();
    s_SphTriProj.resize(2 * s_SphCol.size());
    s_SphColLight.clear();
    s_SphColLight.resize(s_SphCol.size());
}

void ImOrient::CreateArrow()
{
    const int   SUBDIV = 15;
    const float CYL_RADIUS = 0.08f;
    const float CONE_RADIUS = 0.16f;
    const float CONE_LENGTH = 0.25f;
    const float ARROW_BGN = -1.1f;
    const float ARROW_END = 1.15f;
    int i;
    for (i = 0; i < 4; ++i)
    {
        s_ArrowTri[i].clear();
        s_ArrowNorm[i].clear();
    }

    float x0, x1, y0, y1, z0, z1, a0, a1, nx, nn;
    for (i = 0; i < SUBDIV; ++i)
    {
        a0 = 2.0f*glm::pi<float>()*(float(i)) / SUBDIV;
        a1 = 2.0f*glm::pi<float>()*(float(i + 1)) / SUBDIV;
        x0 = ARROW_BGN;
        x1 = ARROW_END - CONE_LENGTH;
        y0 = cosf(a0);
        z0 = sinf(a0);
        y1 = cosf(a1);
        z1 = sinf(a1);
        s_ArrowTri[ARROW_CYL].push_back(x1); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*y0); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*z0);
        s_ArrowTri[ARROW_CYL].push_back(x0); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*y0); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*z0);
        s_ArrowTri[ARROW_CYL].push_back(x0); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*y1); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*z1);
        s_ArrowTri[ARROW_CYL].push_back(x1); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*y0); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*z0);
        s_ArrowTri[ARROW_CYL].push_back(x0); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*y1); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*z1);
        s_ArrowTri[ARROW_CYL].push_back(x1); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*y1); s_ArrowTri[ARROW_CYL].push_back(CYL_RADIUS*z1);
        s_ArrowNorm[ARROW_CYL].push_back(0); s_ArrowNorm[ARROW_CYL].push_back(y0); s_ArrowNorm[ARROW_CYL].push_back(z0);
        s_ArrowNorm[ARROW_CYL].push_back(0); s_ArrowNorm[ARROW_CYL].push_back(y0); s_ArrowNorm[ARROW_CYL].push_back(z0);
        s_ArrowNorm[ARROW_CYL].push_back(0); s_ArrowNorm[ARROW_CYL].push_back(y1); s_ArrowNorm[ARROW_CYL].push_back(z1);
        s_ArrowNorm[ARROW_CYL].push_back(0); s_ArrowNorm[ARROW_CYL].push_back(y0); s_ArrowNorm[ARROW_CYL].push_back(z0);
        s_ArrowNorm[ARROW_CYL].push_back(0); s_ArrowNorm[ARROW_CYL].push_back(y1); s_ArrowNorm[ARROW_CYL].push_back(z1);
        s_ArrowNorm[ARROW_CYL].push_back(0); s_ArrowNorm[ARROW_CYL].push_back(y1); s_ArrowNorm[ARROW_CYL].push_back(z1);
        s_ArrowTri[ARROW_CYL_CAP].push_back(x0); s_ArrowTri[ARROW_CYL_CAP].push_back(0); s_ArrowTri[ARROW_CYL_CAP].push_back(0);
        s_ArrowTri[ARROW_CYL_CAP].push_back(x0); s_ArrowTri[ARROW_CYL_CAP].push_back(CYL_RADIUS*y1); s_ArrowTri[ARROW_CYL_CAP].push_back(CYL_RADIUS*z1);
        s_ArrowTri[ARROW_CYL_CAP].push_back(x0); s_ArrowTri[ARROW_CYL_CAP].push_back(CYL_RADIUS*y0); s_ArrowTri[ARROW_CYL_CAP].push_back(CYL_RADIUS*z0);
        s_ArrowNorm[ARROW_CYL_CAP].push_back(-1); s_ArrowNorm[ARROW_CYL_CAP].push_back(0); s_ArrowNorm[ARROW_CYL_CAP].push_back(0);
        s_ArrowNorm[ARROW_CYL_CAP].push_back(-1); s_ArrowNorm[ARROW_CYL_CAP].push_back(0); s_ArrowNorm[ARROW_CYL_CAP].push_back(0);
        s_ArrowNorm[ARROW_CYL_CAP].push_back(-1); s_ArrowNorm[ARROW_CYL_CAP].push_back(0); s_ArrowNorm[ARROW_CYL_CAP].push_back(0);
        x0 = ARROW_END - CONE_LENGTH;
        x1 = ARROW_END;
        nx = CONE_RADIUS / (x1 - x0);
        nn = 1.0f / sqrtf(nx*nx + 1);
        s_ArrowTri[ARROW_CONE].push_back(x1); s_ArrowTri[ARROW_CONE].push_back(0); s_ArrowTri[ARROW_CONE].push_back(0);
        s_ArrowTri[ARROW_CONE].push_back(x0); s_ArrowTri[ARROW_CONE].push_back(CONE_RADIUS*y0); s_ArrowTri[ARROW_CONE].push_back(CONE_RADIUS*z0);
        s_ArrowTri[ARROW_CONE].push_back(x0); s_ArrowTri[ARROW_CONE].push_back(CONE_RADIUS*y1); s_ArrowTri[ARROW_CONE].push_back(CONE_RADIUS*z1);
        s_ArrowTri[ARROW_CONE].push_back(x1); s_ArrowTri[ARROW_CONE].push_back(0); s_ArrowTri[ARROW_CONE].push_back(0);
        s_ArrowTri[ARROW_CONE].push_back(x0); s_ArrowTri[ARROW_CONE].push_back(CONE_RADIUS*y1); s_ArrowTri[ARROW_CONE].push_back(CONE_RADIUS*z1);
        s_ArrowTri[ARROW_CONE].push_back(x1); s_ArrowTri[ARROW_CONE].push_back(0); s_ArrowTri[ARROW_CONE].push_back(0);
        s_ArrowNorm[ARROW_CONE].push_back(nn*nx); s_ArrowNorm[ARROW_CONE].push_back(nn*y0); s_ArrowNorm[ARROW_CONE].push_back(nn*z0);
        s_ArrowNorm[ARROW_CONE].push_back(nn*nx); s_ArrowNorm[ARROW_CONE].push_back(nn*y0); s_ArrowNorm[ARROW_CONE].push_back(nn*z0);
        s_ArrowNorm[ARROW_CONE].push_back(nn*nx); s_ArrowNorm[ARROW_CONE].push_back(nn*y1); s_ArrowNorm[ARROW_CONE].push_back(nn*z1);
        s_ArrowNorm[ARROW_CONE].push_back(nn*nx); s_ArrowNorm[ARROW_CONE].push_back(nn*y0); s_ArrowNorm[ARROW_CONE].push_back(nn*z0);
        s_ArrowNorm[ARROW_CONE].push_back(nn*nx); s_ArrowNorm[ARROW_CONE].push_back(nn*y1); s_ArrowNorm[ARROW_CONE].push_back(nn*z1);
        s_ArrowNorm[ARROW_CONE].push_back(nn*nx); s_ArrowNorm[ARROW_CONE].push_back(nn*y1); s_ArrowNorm[ARROW_CONE].push_back(nn*z1);
        s_ArrowTri[ARROW_CONE_CAP].push_back(x0); s_ArrowTri[ARROW_CONE_CAP].push_back(0); s_ArrowTri[ARROW_CONE_CAP].push_back(0);
        s_ArrowTri[ARROW_CONE_CAP].push_back(x0); s_ArrowTri[ARROW_CONE_CAP].push_back(CONE_RADIUS*y1); s_ArrowTri[ARROW_CONE_CAP].push_back(CONE_RADIUS*z1);
        s_ArrowTri[ARROW_CONE_CAP].push_back(x0); s_ArrowTri[ARROW_CONE_CAP].push_back(CONE_RADIUS*y0); s_ArrowTri[ARROW_CONE_CAP].push_back(CONE_RADIUS*z0);
        s_ArrowNorm[ARROW_CONE_CAP].push_back(-1); s_ArrowNorm[ARROW_CONE_CAP].push_back(0); s_ArrowNorm[ARROW_CONE_CAP].push_back(0);
        s_ArrowNorm[ARROW_CONE_CAP].push_back(-1); s_ArrowNorm[ARROW_CONE_CAP].push_back(0); s_ArrowNorm[ARROW_CONE_CAP].push_back(0);
        s_ArrowNorm[ARROW_CONE_CAP].push_back(-1); s_ArrowNorm[ARROW_CONE_CAP].push_back(0); s_ArrowNorm[ARROW_CONE_CAP].push_back(0);
    }

    for (i = 0; i < 4; ++i)
    {
        s_ArrowTriProj[i].clear();
        s_ArrowTriProj[i].resize(2 * (s_ArrowTri[i].size() / 3));
        s_ArrowColLight[i].clear();
        s_ArrowColLight[i].resize(s_ArrowTri[i].size() / 3);
    }
}

void ImOrient::Permute(float *outX, float *outY, float *outZ, float x, float y, float z)
{
    float px = x, py = y, pz = z;
    *outX = m_Permute[0][0] * px + m_Permute[1][0] * py + m_Permute[2][0] * pz;
    *outY = m_Permute[0][1] * px + m_Permute[1][1] * py + m_Permute[2][1] * pz;
    *outZ = m_Permute[0][2] * px + m_Permute[1][2] * py + m_Permute[2][2] * pz;
}

void ImOrient::PermuteInv(float *outX, float *outY, float *outZ, float x, float y, float z)
{
    float px = x, py = y, pz = z;
    *outX = m_Permute[0][0]*px + m_Permute[0][1]*py + m_Permute[0][2]*pz;
    *outY = m_Permute[1][0]*px + m_Permute[1][1]*py + m_Permute[1][2]*pz;
    *outZ = m_Permute[2][0]*px + m_Permute[2][1]*py + m_Permute[2][2]*pz;
}

void ImOrient::Permute(double *outX, double *outY, double *outZ, double x, double y, double z)
{
    double px = x, py = y, pz = z;
    *outX = m_Permute[0][0]*px + m_Permute[1][0]*py + m_Permute[2][0]*pz;
    *outY = m_Permute[0][1]*px + m_Permute[1][1]*py + m_Permute[2][1]*pz;
    *outZ = m_Permute[0][2]*px + m_Permute[1][2]*py + m_Permute[2][2]*pz;
}

void ImOrient::PermuteInv(double *outX, double *outY, double *outZ, double x, double y, double z)
{
    double px = x, py = y, pz = z;
    *outX = m_Permute[0][0]*px + m_Permute[0][1]*py + m_Permute[0][2]*pz;
    *outY = m_Permute[1][0]*px + m_Permute[1][1]*py + m_Permute[1][2]*pz;
    *outZ = m_Permute[2][0]*px + m_Permute[2][1]*py + m_Permute[2][2]*pz;
}

void ImOrient::ApplyQuat(float *outX, float *outY, float *outZ, float x, float y, float z, float qx, float qy, float qz, float qs)
{
    float ps = -qx * x - qy * y - qz * z;
    float px = qs * x + qy * z - qz * y;
    float py = qs * y + qz * x - qx * z;
    float pz = qs * z + qx * y - qy * x;
    *outX = -ps * qx + px * qs - py * qz + pz * qy;
    *outY = -ps * qy + py * qs - pz * qx + px * qz;
    *outZ = -ps * qz + pz * qs - px * qy + py * qx;
}

ImOrient::ImOrient()
{
    Qx = Qy = Qz = 0;
    Qs = 1;
    Vx = 1;
    Vy = Vz = 0;
    Angle = 0;
    Dx = Dy = Dz = 0;
    m_AAMode = false; // Axis & angle mode hidden
    m_ShowVal = false;
    m_IsFloat = true;
    m_IsDir = false;
    m_Dir[0] = m_Dir[1] = m_Dir[2] = 0;
    m_DirColor = 0xff00ffff;
    int i, j;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            m_Permute[i][j] = (i == j) ? 1.0f : 0.0f;
    ConvertToAxisAngle();
    m_Highlighted = false;
    m_Rotating = false;
}

/*
void IMGUI_API ImGui_Orient::InitQuat4DCB(void *_ExtValue, void *_ClientData)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(_ExtValue);
    if( ext )
    {
        ext->Qx = ext->Qy = ext->Qz = 0;
        ext->Qs = 1;
        ext->Vx = 1;
        ext->Vy = ext->Vz = 0;
        ext->Angle = 0;
        ext->Dx = ext->Dy = ext->Dz = 0;
        ext->m_AAMode = false; // Axis & angle mode hidden
        ext->m_ShowVal = false;
        ext->m_IsFloat = false;
        ext->m_IsDir = false;
        ext->m_Dir[0] = ext->m_Dir[1] = ext->m_Dir[2] = 0;
        ext->m_DirColor = 0xff00ffff;
        int i, j;
        for(i=0; i<3; ++i)
            for(j=0; j<3; ++j)
                ext->m_Permute[i][j] = (i==j) ? 1.0f : 0.0f;
        ext->m_StructProxy = (CTwMgr::CStructProxy *)_ClientData;
        ext->ConvertToAxisAngle();
        ext->m_Highlighted = false;
        ext->m_Rotating = false;
        if( ext->m_StructProxy!=NULL )
        {
            ext->m_StructProxy->m_CustomDrawCallback = ImGui_Orient::DrawCB;
            ext->m_StructProxy->m_CustomMouseButtonCallback = ImGui_Orient::MouseButtonCB;
            ext->m_StructProxy->m_CustomMouseMotionCallback = ImGui_Orient::MouseMotionCB;
            ext->m_StructProxy->m_CustomMouseLeaveCallback = ImGui_Orient::MouseLeaveCB;
        }
    }
}

void IMGUI_API ImGui_Orient::InitDir3FCB(void *_ExtValue, void *_ClientData)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(_ExtValue);
    if( ext )
    {
        ext->Qx = ext->Qy = ext->Qz = 0;
        ext->Qs = 1;
        ext->Vx = 1;
        ext->Vy = ext->Vz = 0;
        ext->Angle = 0;
        ext->Dx = 1;
        ext->Dy = ext->Dz = 0;
        ext->m_AAMode = false; // Axis & angle mode hidden
        ext->m_ShowVal = true;
        ext->m_IsFloat = true;
        ext->m_IsDir = true;
        ext->m_Dir[0] = ext->m_Dir[1] = ext->m_Dir[2] = 0;
        ext->m_DirColor = 0xff00ffff;
        int i, j;
        for(i=0; i<3; ++i)
            for(j=0; j<3; ++j)
                ext->m_Permute[i][j] = (i==j) ? 1.0f : 0.0f;
        ext->m_StructProxy = (CTwMgr::CStructProxy *)_ClientData;
        ext->ConvertToAxisAngle();
        ext->m_Highlighted = false;
        ext->m_Rotating = false;
        if( ext->m_StructProxy!=NULL )
        {
            ext->m_StructProxy->m_CustomDrawCallback = ImGui_Orient::DrawCB;
            ext->m_StructProxy->m_CustomMouseButtonCallback = ImGui_Orient::MouseButtonCB;
            ext->m_StructProxy->m_CustomMouseMotionCallback = ImGui_Orient::MouseMotionCB;
            ext->m_StructProxy->m_CustomMouseLeaveCallback = ImGui_Orient::MouseLeaveCB;
        }
    }
}

void IMGUI_API ImGui_Orient::InitDir3DCB(void *_ExtValue, void *_ClientData)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(_ExtValue);
    if( ext )
    {
        ext->Qx = ext->Qy = ext->Qz = 0;
        ext->Qs = 1;
        ext->Vx = 1;
        ext->Vy = ext->Vz = 0;
        ext->Angle = 0;
        ext->Dx = 1;
        ext->Dy = ext->Dz = 0;
        ext->m_AAMode = false; // Axis & angle mode hidden
        ext->m_ShowVal = true;
        ext->m_IsFloat = false;
        ext->m_IsDir = true;
        ext->m_Dir[0] = ext->m_Dir[1] = ext->m_Dir[2] = 0;
        ext->m_DirColor = 0xff00ffff;
        int i, j;
        for(i=0; i<3; ++i)
            for(j=0; j<3; ++j)
                ext->m_Permute[i][j] = (i==j) ? 1.0f : 0.0f;
        ext->m_StructProxy = (CTwMgr::CStructProxy *)_ClientData;
        ext->ConvertToAxisAngle();
        ext->m_Highlighted = false;
        ext->m_Rotating = false;
        if( ext->m_StructProxy!=NULL )
        {
            ext->m_StructProxy->m_CustomDrawCallback = ImGui_Orient::DrawCB;
            ext->m_StructProxy->m_CustomMouseButtonCallback = ImGui_Orient::MouseButtonCB;
            ext->m_StructProxy->m_CustomMouseMotionCallback = ImGui_Orient::MouseMotionCB;
            ext->m_StructProxy->m_CustomMouseLeaveCallback = ImGui_Orient::MouseLeaveCB;
        }
    }
}

void IMGUI_API ImGui_Orient::CopyVarFromExtCB(void *_VarValue, const void *_ExtValue, unsigned int _ExtMemberIndex, void *_ClientData)
{
    ImGui_Orient *ext = (ImGui_Orient *)(_ExtValue);
    CTwMgr::CMemberProxy *mProxy = static_cast<CTwMgr::CMemberProxy *>(_ClientData);
    if( _VarValue && ext )
    {
        // Synchronize Quat and AxisAngle
        if( _ExtMemberIndex>=4 && _ExtMemberIndex<=7 )
        {
            ext->ConvertToAxisAngle();
            // show/hide quat values
            if( _ExtMemberIndex==4 && mProxy && mProxy->m_VarParent )
            {
                assert( mProxy->m_VarParent->m_Vars.size()==16 );
                bool visible = ext->m_ShowVal;
                if( ext->m_IsDir )
                {
                    if(    mProxy->m_VarParent->m_Vars[13]->m_Visible != visible
                        || mProxy->m_VarParent->m_Vars[14]->m_Visible != visible
                        || mProxy->m_VarParent->m_Vars[15]->m_Visible != visible )
                    {
                        mProxy->m_VarParent->m_Vars[13]->m_Visible = visible;
                        mProxy->m_VarParent->m_Vars[14]->m_Visible = visible;
                        mProxy->m_VarParent->m_Vars[15]->m_Visible = visible;
                        mProxy->m_Bar->NotUpToDate();
                    }
                }
                else
                {
                    if(    mProxy->m_VarParent->m_Vars[4]->m_Visible != visible
                        || mProxy->m_VarParent->m_Vars[5]->m_Visible != visible
                        || mProxy->m_VarParent->m_Vars[6]->m_Visible != visible
                        || mProxy->m_VarParent->m_Vars[7]->m_Visible != visible )
                    {
                        mProxy->m_VarParent->m_Vars[4]->m_Visible = visible;
                        mProxy->m_VarParent->m_Vars[5]->m_Visible = visible;
                        mProxy->m_VarParent->m_Vars[6]->m_Visible = visible;
                        mProxy->m_VarParent->m_Vars[7]->m_Visible = visible;
                        mProxy->m_Bar->NotUpToDate();
                    }
                }
            }
        }
        else if( _ExtMemberIndex>=8 && _ExtMemberIndex<=11 )
            ext->ConvertFromAxisAngle();
        else if( mProxy && _ExtMemberIndex==12 && mProxy->m_VarParent && !ext->m_IsDir )
        {
            assert( mProxy->m_VarParent->m_Vars.size()==16 );
            bool aa = ext->m_AAMode;
            if(    mProxy->m_VarParent->m_Vars[4]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[5]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[6]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[7]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[8 ]->m_Visible != aa
                || mProxy->m_VarParent->m_Vars[9 ]->m_Visible != aa
                || mProxy->m_VarParent->m_Vars[10]->m_Visible != aa
                || mProxy->m_VarParent->m_Vars[11]->m_Visible != aa )
            {
                mProxy->m_VarParent->m_Vars[4]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[5]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[6]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[7]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[8 ]->m_Visible = aa;
                mProxy->m_VarParent->m_Vars[9 ]->m_Visible = aa;
                mProxy->m_VarParent->m_Vars[10]->m_Visible = aa;
                mProxy->m_VarParent->m_Vars[11]->m_Visible = aa;
                mProxy->m_Bar->NotUpToDate();
            }
            if( static_cast<CTwVarAtom *>(mProxy->m_VarParent->m_Vars[12])->m_ReadOnly )
            {
                static_cast<CTwVarAtom *>(mProxy->m_VarParent->m_Vars[12])->m_ReadOnly = false;
                mProxy->m_Bar->NotUpToDate();
            }
        }

        if( ext->m_IsFloat )
        {
            float *var = static_cast<float *>(_VarValue);
            if( ext->m_IsDir )
            {
                var[0] = (float)ext->Dx;
                var[1] = (float)ext->Dy;
                var[2] = (float)ext->Dz;
            }
            else // quat
            {
                var[0] = (float)ext->Qx;
                var[1] = (float)ext->Qy;
                var[2] = (float)ext->Qz;
                var[3] = (float)ext->Qs;
            }
        }
        else
        {
            double *var = static_cast<double *>(_VarValue);
            if( ext->m_IsDir )
            {
                var[0] = ext->Dx;
                var[1] = ext->Dy;
                var[2] = ext->Dz;
            }
            else // quat
            {
                var[0] = ext->Qx;
                var[1] = ext->Qy;
                var[2] = ext->Qz;
                var[3] = ext->Qs;
            }
        }
    }
}

void IMGUI_API ImGui_Orient::CopyVarToExtCB(const void *_VarValue, void *_ExtValue, unsigned int _ExtMemberIndex, void *_ClientData)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(_ExtValue);
    CTwMgr::CMemberProxy *mProxy = static_cast<CTwMgr::CMemberProxy *>(_ClientData);
    (void)mProxy;
    if( _VarValue && ext )
    {
        if( mProxy && _ExtMemberIndex==12 && mProxy->m_VarParent && !ext->m_IsDir )
        {
            assert( mProxy->m_VarParent->m_Vars.size()==16 );
            bool aa = ext->m_AAMode;
            if(    mProxy->m_VarParent->m_Vars[4]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[5]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[6]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[7]->m_Visible != !aa
                || mProxy->m_VarParent->m_Vars[8 ]->m_Visible != aa
                || mProxy->m_VarParent->m_Vars[9 ]->m_Visible != aa
                || mProxy->m_VarParent->m_Vars[10]->m_Visible != aa
                || mProxy->m_VarParent->m_Vars[11]->m_Visible != aa )
            {
                mProxy->m_VarParent->m_Vars[4]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[5]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[6]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[7]->m_Visible = !aa;
                mProxy->m_VarParent->m_Vars[8 ]->m_Visible = aa;
                mProxy->m_VarParent->m_Vars[9 ]->m_Visible = aa;
                mProxy->m_VarParent->m_Vars[10]->m_Visible = aa;
                mProxy->m_VarParent->m_Vars[11]->m_Visible = aa;
                mProxy->m_Bar->NotUpToDate();
            }
            if( static_cast<CTwVarAtom *>(mProxy->m_VarParent->m_Vars[12])->m_ReadOnly )
            {
                static_cast<CTwVarAtom *>(mProxy->m_VarParent->m_Vars[12])->m_ReadOnly = false;
                mProxy->m_Bar->NotUpToDate();
            }
        }
        else if( mProxy && _ExtMemberIndex==4 && mProxy->m_VarParent )
        {
            assert( mProxy->m_VarParent->m_Vars.size()==16 );
            bool visible = ext->m_ShowVal;
            if( ext->m_IsDir )
            {
                if(    mProxy->m_VarParent->m_Vars[13]->m_Visible != visible
                    || mProxy->m_VarParent->m_Vars[14]->m_Visible != visible
                    || mProxy->m_VarParent->m_Vars[15]->m_Visible != visible )
                {
                    mProxy->m_VarParent->m_Vars[13]->m_Visible = visible;
                    mProxy->m_VarParent->m_Vars[14]->m_Visible = visible;
                    mProxy->m_VarParent->m_Vars[15]->m_Visible = visible;
                    mProxy->m_Bar->NotUpToDate();
                }
            }
            else
            {
                if(    mProxy->m_VarParent->m_Vars[4]->m_Visible != visible
                    || mProxy->m_VarParent->m_Vars[5]->m_Visible != visible
                    || mProxy->m_VarParent->m_Vars[6]->m_Visible != visible
                    || mProxy->m_VarParent->m_Vars[7]->m_Visible != visible )
                {
                    mProxy->m_VarParent->m_Vars[4]->m_Visible = visible;
                    mProxy->m_VarParent->m_Vars[5]->m_Visible = visible;
                    mProxy->m_VarParent->m_Vars[6]->m_Visible = visible;
                    mProxy->m_VarParent->m_Vars[7]->m_Visible = visible;
                    mProxy->m_Bar->NotUpToDate();
                }
            }
        }

        if( ext->m_IsFloat )
        {
            const float *var = static_cast<const float *>(_VarValue);
            if( ext->m_IsDir )
            {
                ext->Dx = var[0];
                ext->Dy = var[1];
                ext->Dz = var[2];
                QuatFromDir(&ext->Qx, &ext->Qy, &ext->Qz, &ext->Qs, var[0], var[1], var[2]);
            }
            else
            {
                ext->Qx = var[0];
                ext->Qy = var[1];
                ext->Qz = var[2];
                ext->Qs = var[3];
            }

        }
        else
        {
            const double *var = static_cast<const double *>(_VarValue);
            if( ext->m_IsDir )
            {
                ext->Dx = var[0];
                ext->Dy = var[1];
                ext->Dz = var[2];
                QuatFromDir(&ext->Qx, &ext->Qy, &ext->Qz, &ext->Qs, var[0], var[1], var[2]);
            }
            else
            {
                ext->Qx = var[0];
                ext->Qy = var[1];
                ext->Qz = var[2];
                ext->Qs = var[3];
            }
        }
        ext->ConvertToAxisAngle();
    }
}
*/
void ImOrient::ConvertToAxisAngle()
{
    if (fabs(Qs) > (1.0 + FLT_EPSILON))
    {
        //Vx = Vy = Vz = 0; // no, keep the previous value
        Angle = 0;
    }
    else
    {
        double a;
        if (Qs >= 1.0f)
            a = 0; // and keep V
        else if (Qs <= -1.0f)
            a = M_PI; // and keep V
        else if (fabs(Qx*Qx + Qy*Qy + Qz*Qz + Qs*Qs) < (FLT_EPSILON * FLT_EPSILON))
            a = 0;
        else
        {
            a = acos(Qs);
            if (a*Angle < 0) // Preserve the sign of Angle
                a = -a;
            double f = 1.0f / sin(a);
            Vx = Qx * f;
            Vy = Qy * f;
            Vz = Qz * f;
        }
        Angle = 2.0*a;
    }

    Angle = RadToDeg(Angle);

    if (fabs(Angle) < FLT_EPSILON && fabs(Vx*Vx + Vy*Vy + Vz*Vz) < FLT_EPSILON * FLT_EPSILON)
        Vx = 1.0e-7;    // all components cannot be null
}

void ImOrient::ConvertFromAxisAngle()
{
    double n = Vx*Vx + Vy*Vy + Vz*Vz;
    if (fabs(n) > (FLT_EPSILON * FLT_EPSILON))
    {
        double f = 0.5*DegToRad(Angle);
        Qs = cos(f);
        f = sin(f);

        Qx = Vx * f;
        Qy = Vy * f;
        Qz = Vz * f;
    }
    else
    {
        Qs = 1.0;
        Qx = Qy = Qz = 0.0;
    }
}

/*
void ImGui_Orient::CopyToVar()
{
    if( m_StructProxy!=NULL )
    {
        if( m_StructProxy->m_StructSetCallback!=NULL )
        {
            if( m_IsFloat )
            {
                if( m_IsDir )
                {
                    float d[] = {1, 0, 0};
                    ApplyQuat(d+0, d+1, d+2, 1, 0, 0, (float)Qx, (float)Qy, (float)Qz, (float)Qs);
                    float l = (float)sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
                    d[0] *= l; d[1] *= l; d[2] *= l;
                    Dx = d[0]; Dy = d[1]; Dz = d[2]; // update also Dx,Dy,Dz
                    m_StructProxy->m_StructSetCallback(d, m_StructProxy->m_StructClientData);
                }
                else
                {
                    float q[] = { (float)Qx, (float)Qy, (float)Qz, (float)Qs };
                    m_StructProxy->m_StructSetCallback(q, m_StructProxy->m_StructClientData);
                }
            }
            else
            {
                if( m_IsDir )
                {
                    float d[] = {1, 0, 0};
                    ApplyQuat(d+0, d+1, d+2, 1, 0, 0, (float)Qx, (float)Qy, (float)Qz, (float)Qs);
                    double l = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
                    double dd[] = {l*d[0], l*d[1], l*d[2]};
                    Dx = dd[0]; Dy = dd[1]; Dz = dd[2]; // update also Dx,Dy,Dz
                    m_StructProxy->m_StructSetCallback(dd, m_StructProxy->m_StructClientData);
                }
                else
                {
                    double q[] = { Qx, Qy, Qz, Qs };
                    m_StructProxy->m_StructSetCallback(q, m_StructProxy->m_StructClientData);
                }
            }
        }
        else if( m_StructProxy->m_StructData!=NULL )
        {
            if( m_IsFloat )
            {
                if( m_IsDir )
                {
                    float *d = static_cast<float *>(m_StructProxy->m_StructData);
                    ApplyQuat(d+0, d+1, d+2, 1, 0, 0, (float)Qx, (float)Qy, (float)Qz, (float)Qs);
                    float l = (float)sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
                    d[0] *= l; d[1] *= l; d[2] *= l;
                    Dx = d[0]; Dy = d[1]; Dz = d[2]; // update also Dx,Dy,Dz
                }
                else
                {
                    float *q = static_cast<float *>(m_StructProxy->m_StructData);
                    q[0] = (float)Qx; q[1] = (float)Qy; q[2] = (float)Qz; q[3] = (float)Qs;
                }
            }
            else
            {
                if( m_IsDir )
                {
                    double *dd = static_cast<double *>(m_StructProxy->m_StructData);
                    float d[] = {1, 0, 0};
                    ApplyQuat(d+0, d+1, d+2, 1, 0, 0, (float)Qx, (float)Qy, (float)Qz, (float)Qs);
                    double l = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
                    dd[0] = l*d[0]; dd[1] = l*d[1]; dd[2] = l*d[2];
                    Dx = dd[0]; Dy = dd[1]; Dz = dd[2]; // update also Dx,Dy,Dz
                }
                else
                {
                    double *q = static_cast<double *>(m_StructProxy->m_StructData);
                    q[0] = Qx; q[1] = Qy; q[2] = Qz; q[3] = Qs;
                }
            }
        }
    }
}
*/
/*
bool ImGui_Orient::MouseMotionCB(int mouseX, int mouseY, int w, int h, void *structExtValue, void *clientData, TwBar *bar, CTwVarGroup *varGrp)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(structExtValue);
    if( ext==NULL )
        return false;
    (void)clientData, (void)varGrp;

    if( mouseX>0 && mouseX<w && mouseY>0 && mouseY<h )
        ext->m_Highlighted = true;

    if( ext->m_Rotating )
    {
        double x = QuatIX(mouseX, w, h);
        double y = QuatIY(mouseY, w, h);
        double z = 1;
        double px, py, pz, ox, oy, oz;
        ext->PermuteInv(&px, &py, &pz, x, y, z);
        ext->PermuteInv(&ox, &oy, &oz, ext->m_OrigX, ext->m_OrigY, 1);
        double n0 = sqrt(ox*ox + oy*oy + oz*oz);
        double n1 = sqrt(px*px + py*py + pz*pz);
        if( n0>DOUBLE_EPS && n1>DOUBLE_EPS )
        {
            double v0[] = { ox/n0, oy/n0, oz/n0 };
            double v1[] = { px/n1, py/n1, pz/n1 };
            double axis[3];
            Vec3Cross(axis, v0, v1);
            double sa = sqrt(Vec3Dot(axis, axis));
            double ca = Vec3Dot(v0, v1);
            double angle = atan2(sa, ca);
            if( x*x+y*y>1.0 )
                angle *= 1.0 + 0.2f*(sqrt(x*x+y*y)-1.0);
            double qrot[4], qres[4], qorig[4];
            QuatFromAxisAngle(qrot, axis, angle);
            double nqorig = sqrt(ext->m_OrigQuat[0]*ext->m_OrigQuat[0]+ext->m_OrigQuat[1]*ext->m_OrigQuat[1]+ext->m_OrigQuat[2]*ext->m_OrigQuat[2]+ext->m_OrigQuat[3]*ext->m_OrigQuat[3]);
            if( fabs(nqorig)>DOUBLE_EPS_SQ )
            {
                qorig[0] = ext->m_OrigQuat[0]/nqorig;
                qorig[1] = ext->m_OrigQuat[1]/nqorig;
                qorig[2] = ext->m_OrigQuat[2]/nqorig;
                qorig[3] = ext->m_OrigQuat[3]/nqorig;
                QuatMult(qres, qrot, qorig);
                ext->Qx = qres[0];
                ext->Qy = qres[1];
                ext->Qz = qres[2];
                ext->Qs = qres[3];
            }
            else
            {
                ext->Qx = qrot[0];
                ext->Qy = qrot[1];
                ext->Qz = qrot[2];
                ext->Qs = qrot[3];
            }
            ext->CopyToVar();
            if( bar!=NULL )
                bar->NotUpToDate();

            ext->m_PrevX = x;
            ext->m_PrevY = y;
        }
    }

    return true;
}

bool ImGui_Orient::MouseButtonCB(TwMouseButtonID button, bool pressed, int mouseX, int mouseY, int w, int h, void *structExtValue, void *clientData, TwBar *bar, CTwVarGroup *varGrp)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(structExtValue);
    if( ext==NULL )
        return false;
    (void)clientData; (void)bar, (void)varGrp;

    if( button==TW_MOUSE_LEFT )
    {
        if( pressed )
        {
            ext->m_OrigQuat[0] = ext->Qx;
            ext->m_OrigQuat[1] = ext->Qy;
            ext->m_OrigQuat[2] = ext->Qz;
            ext->m_OrigQuat[3] = ext->Qs;
            ext->m_OrigX = QuatIX(mouseX, w, h);
            ext->m_OrigY = QuatIY(mouseY, w, h);
            ext->m_PrevX = ext->m_OrigX;
            ext->m_PrevY = ext->m_OrigY;
            ext->m_Rotating = true;
        }
        else
            ext->m_Rotating = false;
    }

    //printf("Click %x\n", structExtValue);
    return true;
}

void ImGui_Orient::MouseLeaveCB(void *structExtValue, void *clientData, TwBar *bar)
{
    ImGui_Orient *ext = static_cast<ImGui_Orient *>(structExtValue);
    if( ext==NULL )
        return;
    (void)clientData; (void)bar;

    //printf("Leave %x\n", structExtValue);
    ext->m_Highlighted = false;
    ext->m_Rotating = false;
}
*/

ImU32 ImOrient::ColorBlend(ImU32 _Color1, ImU32 _Color2, float sigma)
{
    ImColor color1(_Color1);
    ImColor color2(_Color2);
    float invSigma = 1.0f - sigma;

    color1 = ImColor((color1.Value.x * invSigma) + (color2.Value.x * sigma),
        (color1.Value.y * invSigma) + (color2.Value.y * sigma),
        (color1.Value.z * invSigma) + (color2.Value.z * sigma),
        (color1.Value.w * invSigma) + (color2.Value.w * sigma));

    return color1;
}

void ImOrient::QuatMult(double *out, const double *q1, const double *q2)
{
    out[0] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];
    out[1] = q1[3] * q2[1] + q1[1] * q2[3] + q1[2] * q2[0] - q1[0] * q2[2];
    out[2] = q1[3] * q2[2] + q1[2] * q2[3] + q1[0] * q2[1] - q1[1] * q2[0];
    out[3] = q1[3] * q2[3] - (q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2]);
}

void ImOrient::QuatFromAxisAngle(double *out, const double *axis, double angle)
{
    double n = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    if (fabs(n) > DBL_EPSILON)
    {
        double f = 0.5*angle;
        out[3] = cos(f);
        f = sin(f) / sqrt(n);
        out[0] = axis[0] * f;
        out[1] = axis[1] * f;
        out[2] = axis[2] * f;
    }
    else
    {
        out[3] = 1.0;
        out[0] = out[1] = out[2] = 0.0;
    }
}

void ImOrient::QuatFromDir(double *outQx, double *outQy, double *outQz, double *outQs, double dx, double dy, double dz)
{
    // compute a quaternion that rotates (1,0,0) to (dx,dy,dz)
    double dn = sqrt(dx*dx + dy*dy + dz*dz);
    if (dn < DBL_EPSILON * DBL_EPSILON)
    {
        *outQx = *outQy = *outQz = 0;
        *outQs = 1;
    }
    else
    {
        double rotAxis[3] = { 0, -dz, dy };
        if (rotAxis[0] * rotAxis[0] + rotAxis[1] * rotAxis[1] + rotAxis[2] * rotAxis[2] < (DBL_EPSILON * DBL_EPSILON))
        {
            rotAxis[0] = rotAxis[1] = 0;
            rotAxis[2] = 1;
        }
        double rotAngle = acos(dx / dn);
        double rotQuat[4];
        QuatFromAxisAngle(rotQuat, rotAxis, rotAngle);
        *outQx = rotQuat[0];
        *outQy = rotQuat[1];
        *outQz = rotQuat[2];
        *outQs = rotQuat[3];
    }
}

void ImOrient::Vec3Cross(double *out, const double *a, const double *b)
{
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

double ImOrient::Vec3Dot(const double *a, const double *b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void ImOrient::Vec3RotY(float *x, float *y, float *z)
{
    (void)y;
    float tmp = *x;
    *x = -*z;
    *z = tmp;
}

void ImOrient::Vec3RotZ(float *x, float *y, float *z)
{
    (void)z;
    float tmp = *x;
    *x = -*y;
    *y = tmp;
}

