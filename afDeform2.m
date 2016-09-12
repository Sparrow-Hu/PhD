%this is to show the results after optimization and plot them
% for 2 patches
function afDeform2(X)
global  var_ks_g1_u var_ks_g1_v para_cstr_patch1 para_cstr_patch2 tLoop1 tLoop2   num_Lp1_cstr    var_ks_g2_u var_ks_g2_v num_Lp2_cstr  ...
    start_p1 u_end1 v_end1 end_p1 num_u1 num_v1   start_p2 u_end2 v_end2 end_p2 num_u2 num_v2  var_cp_g1 var_ws_g1  var_cp_p1 var_ws_p1 var_ks_p1 var_cp_g2 ...
    var_ws_g2 var_cp_p2 var_ws_p2 var_ks_p2 l_x0 l_x1  surf1 surf2 ;%this flag is to control the cases of 3D/2D/3_2D
%% define the deformed surfaces
x0=[];
x1=[];
j=1;
for i=1:l_x0
    x0(j) = X(i);
    j=j+1;
end
j=1;
for i=l_x0+1:l_x0+l_x1
    x1(j) = X(i);
    j=j+1;
end
% j=1;
% for i=l_x0+l_x1+1:l_x0+l_x1+l_x2
%     x2(j) = X(i);
%     j=j+1;
% end
% j=1;
% for i=l_x0+l_x1+l_x2+1:l_x0+l_x1+l_x2+l_x3
%     x3(j) = X(i);
%     j=j+1;
% end
if isempty(x0)
    srfc_1 = surf1;%if there is no variables for deformation surface, use the original surface
else
    [x_g1,y_g1,z_g1,h_g1,k_g1_u,k_g1_v,x_p1,y_p1,h_p1,k_p1]=paraSplit(x0,var_cp_g1,var_ws_g1,var_ks_g1_u,var_ks_g1_v,var_cp_p1,var_ws_p1,var_ks_p1);
    %construct the deformed points of patch 1
    srfc_1= deformSurf(start_p1,u_end1,v_end1,end_p1,num_u1,num_v1,x_g1,y_g1,z_g1,var_cp_g1,h_g1,var_ws_g1,k_g1_u,k_g1_v,var_ks_g1_u,var_ks_g1_v,surf1);%construct the deformed patch1
end
if isempty(x1)
    srfc_2 = surf2;%if there is no variables for deformation surface, use the original surface
else
    [x_g2,y_g2,z_g2,h_g2,k_g2_u,k_g2_v,x_p2,y_p2,h_p2,k_p2]=paraSplit(x1,var_cp_g2,var_ws_g2,var_ks_g2_u,var_ks_g2_v,var_cp_p2,var_ws_p2,var_ks_p2);
    srfc_2= deformSurf(start_p2,u_end2,v_end2,end_p2,num_u2,num_v2,x_g2,y_g2,z_g2,var_cp_g2,h_g2,var_ws_g2,k_g2_u,k_g2_v,var_ks_g2_u,var_ks_g2_v,surf2);% construct the deformed patch2
end
% if isempty(x2)
%     srfc_3 = surf3;%if there is no variables for deformation surface, use the original surface
% else
%     [x_g3,y_g3,z_g3,h_g3,k_g3_u,k_g3_v,x_p3,y_p3,h_p3,k_p3]=paraSplit(x2,var_cp_g3,var_ws_g3,var_ks_g3_u,var_ks_g3_v,var_cp_p3,var_ws_p3,var_ks_p3);
%     %construct the deformed points of patch 1
%     srfc_3= deformSurf(start_p3,u_end3,v_end3,end_p3,num_u3,num_v3,x_g3,y_g3,z_g3,var_cp_g3,h_g3,var_ws_g3,k_g3_u,k_g3_v,var_ks_g3_u,var_ks_g3_v,surf3);%construct the deformed patch1
% end
% if isempty(x3)
%     srfc_4 = surf4;%if there is no variables for deformation surface, use the original surface
% else
%     [x_g4,y_g4,z_g4,h_g4,k_g4_u,k_g4_v,x_p4,y_p4,h_p4,k_p4]=paraSplit(x3,var_cp_g4,var_ws_g4,var_ks_g4_u,var_ks_g4_v,var_cp_p4,var_ws_p4,var_ks_p4);
%     srfc_4= deformSurf(start_p4,u_end4,v_end4,end_p4,num_u4,num_v4,x_g4,y_g4,z_g4,var_cp_g4,h_g4,var_ws_g4,k_g4_u,k_g4_v,var_ks_g4_u,var_ks_g4_v,surf4);% construct the deformed patch2
% end
% define the deformed trimming curve in parameric space
% for trimming curve of patch 1
%%
if ~isempty(tLoop1)
    trimCrv_p1 = TrimCurve(tLoop1{num_Lp1_cstr});% initialize trimming curve
end
if ~isempty(var_cp_p1)
    for i=1:length(var_cp_p1)% this is to set the control points of trimming curve of patch1 as variables
        cp = [var_cp_p1(i),x_p1(i),y_p1(i)];
        trimCrv_p1.CtrPts(cp);
    end
end
if ~isempty(var_ws_p1)
    for i=1:length(var_ws_p1)% this is to set the weights of trimming curve of patch1 as variables
        ws = [var_ws_p1(i),h_p1(i)];
        trimCrv_p1.Weights(ws);
    end
end
if ~isempty(var_ks_p1)
    for i=1:length(var_ks_p1)% this is to set the knots of trimming curve of patch1 as variables
        ks = [var_ks_p1(i),k_p1(i)];
        trimCrv_p1.Knots(ks);
    end
end
% for trimming curve of patch 2
if ~isempty(tLoop2)
    trimCrv_p2 = TrimCurve(tLoop2{num_Lp2_cstr});%initialize the trimming curve
end
if ~isempty(var_cp_p2)
    for i=1:length(var_cp_p2)% this is to set the control points of trimming curve of patch2 as variables
        cp = [var_cp_p2(i),x_p2(i),y_p2(i)];
        trimCrv_p2.CtrPts(cp);
    end
end
if ~isempty(var_ws_p2)
    for i=1:length(var_ws_p2)% this is to set the weights of trimming curve of patch2 as variables
        ws = [var_ws_p2(i),h_p2(i)];
        trimCrv_p2.Weights(ws);
    end
end
if ~isempty(var_ks_p2)
    for i=1:length(var_ks_p2)% this is to set the knots of trimming curve of patch1 as variables
        ks = [var_ks_p2(i),k_p2(i)];
        trimCrv_p2.Knots(ks);
    end
end
%%
%reconstruct the trimming loops
if ~isempty(tLoop1)
    [m1,n1]=size(tLoop1{1}.coefs);% for trimming loop1
    tLoop1{1}.coefs(:,n1) = trimCrv_p1.ini_value.coefs(:,1);%maintain the connection between four trimming curves
    [m2,n2]=size(trimCrv_p1.ini_value.coefs);% for trimming loop1
    tLoop1{3}.coefs(:,1) = trimCrv_p1.ini_value.coefs(:,n2);%maintain the connection between four trimming curves
    tLoop1{2} = trimCrv_p1.ini_value;
    plotTrimCrv3D(srfc_1,tLoop1);
hold on;
end
% for trimming loop2
if ~isempty(tLoop2)
    [m3,n3]=size(trimCrv_p2.ini_value.coefs);% for trimming curve 4
    tLoop2{1}.coefs(:,1) = trimCrv_p2.ini_value.coefs(:,n3);
    [m4,n4]=size(tLoop2{3}.coefs);% for trimming curve 4
    tLoop2{3}.coefs(:,n4) = trimCrv_p2.ini_value.coefs(:,1);
    tLoop2{4} = trimCrv_p2.ini_value;
    plotTrimCrv3D(srfc_2,tLoop2);
hold on;
end
%%
if ~isempty(tLoop2) && ~isempty(tLoop1)
    plotNURBS(srfc_1);plotTrimCrv3D(srfc_1,tLoop1);hold on;
    plotNURBS(srfc_2,[0.5,0.5,0]);hold on;plotTrimCrv3D(srfc_2,tLoop2);hold off;
end
if ~isempty(para_cstr_patch1) && ~isempty(para_cstr_patch2)
plotNURBS(srfc_1);hold on;plotParapts(srfc_1,para_cstr_patch1);hold on;
plotNURBS(srfc_2,[0.5,0.5,0]);hold on;plotParapts(srfc_2,para_cstr_patch2);hold off;
% subplot(2,2,3);
% plotNURBS(tLoop1{1});hold on;
% plotNURBS(tLoop1{2});hold on;
% plotNURBS(tLoop1{3});hold on;
% plotNURBS(tLoop1{4});hold off;
% %%plot trimming loop
% subplot(2,2,4);
% plotNURBS(tLoop2{1});hold on;
% plotNURBS(tLoop2{2});hold on;
% plotNURBS(tLoop2{3});hold on;
% plotNURBS(tLoop2{4});hold off;
end