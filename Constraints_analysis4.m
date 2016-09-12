function [ ceq ] = Constraints_analysis4(X)
global    surf3 surf4 var_cp_g3 var_ws_g3  var_cp_p3 var_ws_p3 var_ks_p3 var_cp_g4 var_ws_g4  var_cp_p4 var_ws_p4 var_ks_p4 var_ks_g1_u var_ks_g1_v var_ks_g4_u var_ks_g4_v...
    tLoop1 tLoop2  trimStatus_cstr patchStatus_cstr num_Lp1_cstr para_cstr_trim1 para_cstr_trim2  para_cstr_patch3 para_cstr_patch4 var_ks_g2_u var_ks_g2_v var_ks_g3_u var_ks_g3_v...
    num_Lp2_cstr pts_cstr_trim_p1 pts_cstr_trim_p2 para_cstr_patch1 para_cstr_patch2  num_cstr_trim num_cstr_patch cstrMatrix_trim cstrMatrix_patch  ...
    start_p1 u_end1 v_end1 end_p1 num_u1 num_v1  start_p3 u_end3 v_end3 end_p3 num_u3 num_v3 start_p2 u_end2 v_end2 end_p2 num_u2 num_v2 start_p4 u_end4 v_end4 end_p4 num_u4 num_v4 var_cp_g1 var_ws_g1 var_ks_g1 var_cp_p1 var_ws_p1 var_ks_p1 var_cp_g2 ...
    var_ws_g2 var_cp_p2 var_ws_p2 var_ks_p2 l_x0 l_x1 l_x2 l_x3 surf1 surf2 ;
%% define the deformed surfaces
ceq=[];
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
j=1;
for i=l_x0+l_x1+1:l_x0+l_x1+l_x2
    x2(j) = X(i);
    j=j+1;
end
j=1;
for i=l_x0+l_x1+l_x2+1:l_x0+l_x1+l_x2+l_x3
    x3(j) = X(i);
    j=j+1;
end
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
if isempty(x2)
    srfc_3 = surf3;%if there is no variables for deformation surface, use the original surface
else
    [x_g3,y_g3,z_g3,h_g3,k_g3_u,k_g3_v,x_p3,y_p3,h_p3,k_p3]=paraSplit(x2,var_cp_g3,var_ws_g3,var_ks_g3_u,var_ks_g3_v,var_cp_p3,var_ws_p3,var_ks_p3);
    %construct the deformed points of patch 1
    srfc_3= deformSurf(start_p3,u_end3,v_end3,end_p3,num_u3,num_v3,x_g3,y_g3,z_g3,var_cp_g3,h_g3,var_ws_g3,k_g3_u,k_g3_v,var_ks_g3_u,var_ks_g3_v,surf3);%construct the deformed patch1
end
if isempty(x3)
    srfc_4 = surf4;%if there is no variables for deformation surface, use the original surface
else
    [x_g4,y_g4,z_g4,h_g4,k_g4_u,k_g4_v,x_p4,y_p4,h_p4,k_p4]=paraSplit(x3,var_cp_g4,var_ws_g4,var_ks_g4_u,var_ks_g4_v,var_cp_p4,var_ws_p4,var_ks_p4);
    srfc_4= deformSurf(start_p4,u_end4,v_end4,end_p4,num_u4,num_v4,x_g4,y_g4,z_g4,var_cp_g4,h_g4,var_ws_g4,k_g4_u,k_g4_v,var_ks_g4_u,var_ks_g4_v,surf4);% construct the deformed patch2
end

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
%%  below is to use deformed surfaces to calculate constraints
%
geops_cstr_trim1 = [];geops_cstr_trim2 = [];geops_cstr_patch1 = [];geops_cstr_patch2 = [];geops_cstr_patch3 = [];geops_cstr_patch4 = [];
geops_cstr_trim1_1 = [];geops_cstr_trim2_1 = [];geops_cstr_patch1_1 = [];geops_cstr_patch2_1 = [];geops_cstr_patch3_1 = [];geops_cstr_patch4_1 = [];
if isTrue(trimStatus_cstr) && ~isTrue(patchStatus_cstr)
    for i=1:num_cstr_trim
        if isempty(para_cstr_trim1) && ~isempty(pts_cstr_trim_p1)
            if isym(trimCeval(trimCrv_p1.ini_value,pts_cstr_trim_p1(i)))
                para_cstr_trim1 = sym(para_cstr_trim1);%symbolic representation
            end
            para_cstr_trim1(i,:) = trimCeval(trimCrv_p1.ini_value,pts_cstr_trim_p1(i)); % store the constraint points in parametric space---trimmming loop1
        end
        pct1=para_cstr_trim1(i,:);
        if isym(nrbeval(srfc_1,{pct1(1),pct1(2)}))
            geops_cstr_trim1 = sym(geops_cstr_trim1);%symbolic representation
        end
        geops_cstr_trim1(i,:) = nrbeval(srfc_1,{pct1(1),pct1(2)});
        if isempty(para_cstr_trim2) && ~isempty(pts_cstr_trim_p2)
            if isym(trimCeval(trimCrv_p2.ini_value,pts_cstr_trim_p2(i)))
                para_cstr_trim2=sym(para_cstr_trim2);%symbolic representation
            end% store the constraint points in geometric space
            para_cstr_trim2(i,:) = trimCeval(trimCrv_p2.ini_value,pts_cstr_trim_p2(i)); % store the constraint points in parametric space---trimmming loop2
        end
        pct2=para_cstr_trim2(i,:);
        if isym(nrbeval(srfc_2,{pct2(1),pct2(2)}))
            geops_cstr_trim2=sym(geops_cstr_trim2);%symbolic representation
        end
        geops_cstr_trim2(i,:) = nrbeval(srfc_2,{pct2(1),pct2(2)});            % store the constraint points in geometric space
    end
end
if isTrue(patchStatus_cstr) && ~isTrue(trimStatus_cstr)
    for i=1:num_cstr_patch
        pcp1 = para_cstr_patch1(i,:);% parameters on side 1
        pcp2 = para_cstr_patch2(i,:);% parameters on side 2
        pcp3 = para_cstr_patch3(i,:);% parameters on side 3
        pcp4 = para_cstr_patch4(i,:);% parameters on side 4
        if isym(nrbeval(srfc_1,{pcp1(1),pcp1(2)}))
            geops_cstr_patch1=sym(geops_cstr_patch1);
        end
        if isym(nrbeval(srfc_1,{pcp3(1),pcp3(2)}))
            geops_cstr_patch1_1=sym(geops_cstr_patch1_1);
        end
        geops_cstr_patch1(i,:) = nrbeval(srfc_1,{pcp1(1),pcp1(2)});% side 1/patch 1
        geops_cstr_patch1_1(i,:) = nrbeval(srfc_1,{pcp3(1),pcp3(2)}); % side 3/patch 1
        if isym(nrbeval(srfc_2,{pcp2(1),pcp2(2)}))
            geops_cstr_patch2=sym(geops_cstr_patch2);
        end
        if isym(nrbeval(srfc_2,{pcp3(1),pcp3(2)}))
            geops_cstr_patch2_1=sym(geops_cstr_patch2_1);
        end
        geops_cstr_patch2(i,:) = nrbeval(srfc_2,{pcp2(1),pcp2(2)}); % side 2/patch 2
        geops_cstr_patch2_1(i,:) = nrbeval(srfc_2,{pcp3(1),pcp3(2)}); % side 3/patch 2
        if isym(nrbeval(srfc_3,{pcp4(1),pcp4(2)}))
            geops_cstr_patch3=sym(geops_cstr_patch3);
        end
        if isym(nrbeval(srfc_3,{pcp1(1),pcp1(2)}))
            geops_cstr_patch3_1=sym(geops_cstr_patch3_1);
        end
        geops_cstr_patch3(i,:) = nrbeval(srfc_3,{pcp4(1),pcp4(2)}); % side 4/patch 3
        geops_cstr_patch3_1(i,:) = nrbeval(srfc_3,{pcp1(1),pcp1(2)}); % side 1/patch 3
        if isym(nrbeval(srfc_4,{pcp4(1),pcp4(2)}))
            geops_cstr_patch4=sym(geops_cstr_patch4);
        end
        if isym(nrbeval(srfc_4,{pcp2(1),pcp2(2)}))
            geops_cstr_patch4_1=sym(geops_cstr_patch4_1);
        end
        geops_cstr_patch4(i,:) = nrbeval(srfc_4,{pcp4(1),pcp4(2)}); % side 4/patch 4
        geops_cstr_patch4_1(i,:) = nrbeval(srfc_4,{pcp2(1),pcp2(2)}); % side 2/patch 4
    end
end
if isTrue(patchStatus_cstr) && isTrue(trimStatus_cstr)
    for i=1:num_cstr_trim
        if isempty(para_cstr_trim1) && ~isempty(pts_cstr_trim_p1)
            if isym(trimCeval(trimCrv_p1.ini_value,pts_cstr_trim_p1(i)))
                para_cstr_trim1 = sym(para_cstr_trim1);%symbolic representation
            end
            para_cstr_trim1(i,:) = trimCeval(trimCrv_p1.ini_value,pts_cstr_trim_p1(i)); % store the constraint points in parametric space---trimmming loop1
        end
        pct1=para_cstr_trim1(i,:);
        if isym(nrbeval(srfc_1,{pct1(1),pct1(2)}))
            geops_cstr_trim1 = sym(geops_cstr_trim1);%symbolic representation
        end
        geops_cstr_trim1(i,:) = nrbeval(srfc_1,{pct1(1),pct1(2)});
        if isempty(para_cstr_trim2) && ~isempty(pts_cstr_trim_p2)
            if isym(trimCeval(trimCrv_p2.ini_value,pts_cstr_trim_p2(i)))
                para_cstr_trim2=sym(para_cstr_trim2);%symbolic representation
            end% store the constraint points in geometric space
            para_cstr_trim2(i,:) = trimCeval(trimCrv_p2.ini_value,pts_cstr_trim_p2(i)); % store the constraint points in parametric space---trimmming loop2
        end
        pct2=para_cstr_trim2(i,:);
        if isym(nrbeval(srfc_2,{pct2(1),pct2(2)}))
            geops_cstr_trim2=sym(geops_cstr_trim2);%symbolic representation
        end
        geops_cstr_trim2(i,:) = nrbeval(srfc_2,{pct2(1),pct2(2)});            % store the constraint points in geometric space
    end
    for i=1:num_cstr_patch
        pcp1 = para_cstr_patch1(i,:);
        if isym(nrbeval(srfc_1,{pcp1(1),pcp1(2)}))
            geops_cstr_patch1=sym(geops_cstr_patch1);
        end
        geops_cstr_patch1(i,:) = nrbeval(srfc_1,{pcp1(1),pcp1(2)}); % store the constraints on patches in geometric space
        pcp2 = para_cstr_patch2(i,:);
        if isym(nrbeval(srfc_2,{pcp2(1),pcp2(2)}))
            geops_cstr_patch2=sym(geops_cstr_patch2);
        end
        geops_cstr_patch2(i,:) = nrbeval(srfc_2,{pcp2(1),pcp2(2)});           % store the constraints on patches in geometric space
    end
end
   
temp=0;
%compute the constraints between pair of constraint points
if isTrue(patchStatus_cstr) && ~isTrue(trimStatus_cstr)
    %% constraints between side 1 patch 1 and side 2 patch 2
    for i = 1: num_cstr_patch
        if cstrMatrix_patch(i,1) == 1 %coincide constraints
            if isym(geops_cstr_patch1)  |  isym(geops_cstr_patch2)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_patch1(i,:) - geops_cstr_patch2(i,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,2) == 1  %position constraints
            if isym(norm(geops_cstr_patch1(i,:),geops_cstr_patch2(i,:)))
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_patch1(i,:),geops_cstr_patch2(i,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,3) == 1 %tangent constraints
            pcp1=para_cstr_patch1(i,:);
            [pt1,dpt1] = nrbdeval(srfc_1,nrbderiv(srfc_1),{pcp1(1),pcp1(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pcp2=para_cstr_patch2(i,:);
            [pt2,dpt2] = nrbdeval(srfc_2,nrbderiv(srfc_2),{pcp2(1),pcp2(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_patch(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
    %% constraints between s3p1 and s4p3
    for i = num_cstr_patch+1:num_cstr_patch*2
        if cstrMatrix_patch(i,1) == 1 %coincide constraints
            if isym(geops_cstr_patch1_1)  |  isym(geops_cstr_patch3)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_patch1_1(i-num_cstr_patch,:) - geops_cstr_patch3(i-num_cstr_patch,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,2) == 1  %position constraints
            if isym(norm(geops_cstr_patch1_1(i-num_cstr_patch,:),geops_cstr_patch3(i-num_cstr_patch,:)))
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_patch1_1(i-num_cstr_patch,:),geops_cstr_patch3(i-num_cstr_patch,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,3) == 1 %tangent constraints
            pcp3=para_cstr_patch3(i-num_cstr_patch,:);
            [pt1,dpt1] = nrbdeval(srfc_1,nrbderiv(srfc_1),{pcp3(1),pcp3(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pcp4=para_cstr_patch4(i-num_cstr_patch,:);
            [pt2,dpt2] = nrbdeval(srfc_3,nrbderiv(srfc_3),{pcp4(1),pcp4(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_patch(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
    %% constraints between s1p3 and s2p4
    for i = num_cstr_patch*2+1:num_cstr_patch*3
        if cstrMatrix_patch(i,1) == 1 %coincide constraints
            if isym(geops_cstr_patch3_1)  |  isym(geops_cstr_patch4_1)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_patch3_1(i-num_cstr_patch*2,:) - geops_cstr_patch4_1(i-num_cstr_patch*2,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,2) == 1  %position constraints
            if isym(norm(geops_cstr_patch3_1(i-num_cstr_patch*2,:),geops_cstr_patch4_1(i-num_cstr_patch*2,:)))
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_patch3_1(i,:),geops_cstr_patch4_1(i,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,3) == 1 %tangent constraints
            pcp1=para_cstr_patch1(i-num_cstr_patch*2,:);
            [pt1,dpt1] = nrbdeval(srfc_3,nrbderiv(srfc_3),{pcp1(1),pcp1(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pcp2=para_cstr_patch2(i-num_cstr_patch*2,:);
            [pt2,dpt2] = nrbdeval(srfc_4,nrbderiv(srfc_4),{pcp2(1),pcp2(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_patch(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
    %% constraints between s3p2 and s4p4
    for i = num_cstr_patch*3+1:num_cstr_patch*4
        if cstrMatrix_patch(i,1) == 1 %coincide constraints
            if isym(geops_cstr_patch4)  |  isym(geops_cstr_patch2_1)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_patch4(i-num_cstr_patch*3,:) - geops_cstr_patch2_1(i-num_cstr_patch*3,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,2) == 1  %position constraints
            if isym(norm(geops_cstr_patch4(i-num_cstr_patch*3,:),geops_cstr_patch2_1(i-num_cstr_patch*3,:)))
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_patch4(i-num_cstr_patch*3,:),geops_cstr_patch2_1(i-num_cstr_patch*3,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,3) == 1 %tangent constraints
            pcp4=para_cstr_patch4(i-num_cstr_patch*3,:);
            [pt1,dpt1] = nrbdeval(srfc_4,nrbderiv(srfc_4),{pcp4(1),pcp4(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pcp3=para_cstr_patch3(i-num_cstr_patch*3,:);
            [pt2,dpt2] = nrbdeval(srfc_2,nrbderiv(srfc_2),{pcp3(1),pcp3(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_patch(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
end
if isTrue(trimStatus_cstr) && ~isTrue(patchStatus_cstr)
    for i = 1: num_cstr_trim
        if cstrMatrix_trim(i,1) == 1 %coincide constraints
            if isym(geops_cstr_trim1) | isym(geops_cstr_trim2)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_trim1(i,:) - geops_cstr_trim2(i,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_trim(i,2) == 1  %position constraints
            if isym(geops_cstr_trim1) | isym(geops_cstr_trim2)
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_trim1(i,:),geops_cstr_trim2(i,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_trim(i,3) == 1 %tangent constraints
            pct1=para_cstr_trim1(i,:);
            [pt1,dpt1] = nrbdeval(srfc_1,nrbderiv(srfc_1),{pct1(1),pct1(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pct2=para_cstr_trim2(i,:);
            [pt2,dpt2] = nrbdeval(srfc_2,nrbderiv(srfc_2),{pct2(1),pct2(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_trim(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end  
end
if isTrue(trimStatus_cstr) && ~isTrue(patchStatus_cstr)
    for i = 1: num_cstr_trim
        if cstrMatrix_trim(i,1) == 1 %coincide constraints
            if isym(geops_cstr_trim1) | isym(geops_cstr_trim2)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_trim1(i,:) - geops_cstr_trim2(i,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_trim(i,2) == 1  %position constraints
            if isym(geops_cstr_trim1) | isym(geops_cstr_trim2)
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_trim1(i,:),geops_cstr_trim2(i,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_trim(i,3) == 1 %tangent constraints
            pct1=para_cstr_trim1(i,:);
            [pt1,dpt1] = nrbdeval(srfc_1,nrbderiv(srfc_1),{pct1(1),pct1(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pct2=para_cstr_trim2(i,:);
            [pt2,dpt2] = nrbdeval(srfc_2,nrbderiv(srfc_2),{pct2(1),pct2(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_trim(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
    
end
if isTrue(patchStatus_cstr) && isTrue(trimStatus_cstr)
    for i = 1: num_cstr_trim%%constraints on trimming curves
        if cstrMatrix_trim(i,1) == 1 %coincide constraints
            if isym(geops_cstr_trim1) | isym(geops_cstr_trim2)
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_trim1(i,:) - geops_cstr_trim2(i,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_trim(i,2) == 1  %position constraints
            if isym(geops_cstr_trim1) | isym(geops_cstr_trim2)
                ceq=sym(ceq);
            end
            ceq(temp+1) = norm(geops_cstr_trim1(i,:),geops_cstr_trim2(i,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_trim(i,3) == 1 %tangent constraints
            pct1=para_cstr_trim1(i,:);
            [pt1,dpt1] = nrbdeval(srfc_1,nrbderiv(srfc_1),{pct1(1),pct1(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pct2=para_cstr_trim2(i,:);
            [pt2,dpt2] = nrbdeval(srfc_2,nrbderiv(srfc_2),{pct2(1),pct2(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_trim(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
    for i = 1: num_cstr_patch%%second,constraints points on patches
        if cstrMatrix_patch(i,1) == 1 %
            if isym(geops_cstr_patch1)  |  isym(geops_cstr_patch2)
               ceq=sym(ceq);
            end
            ceq(temp+1:temp+3) = geops_cstr_patch1(i,:) - geops_cstr_patch2(i,:);
            temp =temp+3;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,2) == 1  %position constraints
            if isym(norm(geops_cstr_patch1(i,:),geops_cstr_patch2(i,:)))
              ceq=sym(ceq);
            end  
            ceq(temp+1) = norm(geops_cstr_patch1(i,:),geops_cstr_patch2(i,:));
            temp = temp+1;
        else
            temp = temp;
        end
        if cstrMatrix_patch(i,3) == 1 %tangent constraints
            pcp1=para_cstr_patch1(i,:);
            [pt1,dpt1] = nrbdeval(srfc_1,nrbderiv(srfc_1),{pcp1(1),pcp1(2)});
            norm1=cross(dpt1{1},dpt1{2});
            pcp2=para_cstr_patch2(i,:);
            [pt2,dpt2] = nrbdeval(srfc_2,nrbderiv(srfc_2),{pcp2(1),pcp2(2)});
            norm2 = cross(dpt2{1},dpt2{2});
            if isym(cross(norm1,norm2))
                ceq=sym(ceq);
            end
            ceq(temp+1:temp+3)= cross(norm1,norm2);
            temp=temp+3;
        else
            temp =temp;
        end
        if cstrMatrix_patch(i,4) == 1 %curvature constraints
            disp('this is to define the curvature constraints,but not yet defined!')
        end
        
    end
    
    
end



end

      
       
      

