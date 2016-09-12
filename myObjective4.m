function   [f]= myObjective4(X)
global   var_ks_g1_u var_ks_g1_v  surf3 surf4 var_cp_g3 var_ws_g3  var_cp_p3 var_ws_p3 var_ks_p3 var_cp_g4 var_ws_g4  var_cp_p4 var_ws_p4 var_ks_p4 tLoop1 tLoop2  var_ks_g2_u var_ks_g2_v num_Lp1_cstr var_ks_g3_u var_ks_g3_v obj var_ks_g4_u var_ks_g4_v num_Lp2_cstr ...
    start_p1 u_end1 v_end1 end_p1 num_u1 num_v1  start_p3 u_end3 v_end3 end_p3 num_u3 num_v3 start_p2 u_end2 v_end2 end_p2 num_u2 num_v2 start_p4 u_end4 v_end4 end_p4 num_u4 num_v4 var_cp_g1 var_ws_g1 var_cp_p1 var_ws_p1 var_ks_p1 var_cp_g2 var_ws_g2  var_cp_p2 var_ws_p2 var_ks_p2  l_x0 l_x1 l_x2 l_x3 surf1 surf2  ;
%c=[];
%% below is to define the deformed surfaces
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
    srfo_1 = surf1;%if there is no variables for deformation surface, use the original surface
else
    [x_g1,y_g1,z_g1,h_g1,k_g1_u,k_g1_v,x_p1,y_p1,h_p1,k_p1]=paraSplit(x0,var_cp_g1,var_ws_g1,var_ks_g1_u,var_ks_g1_v,var_cp_p1,var_ws_p1,var_ks_p1);
    %construct the deformed points of patch 1
    srfo_1= deformSurf(start_p1,u_end1,v_end1,end_p1,num_u1,num_v1,x_g1,y_g1,z_g1,var_cp_g1,h_g1,var_ws_g1,k_g1_u,k_g1_v,var_ks_g1_u,var_ks_g1_v,surf1);%construct the deformed patch1
end
if isempty(x1)
    srfo_2 = surf2;%if there is no variables for deformation surface, use the original surface
else
    [x_g2,y_g2,z_g2,h_g2,k_g2_u,k_g2_v,x_p2,y_p2,h_p2,k_p2]=paraSplit(x1,var_cp_g2,var_ws_g2,var_ks_g2_u,var_ks_g2_v,var_cp_p2,var_ws_p2,var_ks_p2);
    srfo_2= deformSurf(start_p2,u_end2,v_end2,end_p2,num_u2,num_v2,x_g2,y_g2,z_g2,var_cp_g2,h_g2,var_ws_g2,k_g2_u,k_g2_v,var_ks_g2_u,var_ks_g2_v,surf2);% construct the deformed patch2
end
if isempty(x2)
    srfo_3 = surf3;%if there is no variables for deformation surface, use the original surface
else
    [x_g3,y_g3,z_g3,h_g3,k_g3_u,k_g3_v,x_p3,y_p3,h_p3,k_p3]=paraSplit(x2,var_cp_g3,var_ws_g3,var_ks_g3_u,var_ks_g3_v,var_cp_p3,var_ws_p3,var_ks_p3);
    %construct the deformed points of patch 1
    srfo_3= deformSurf(start_p3,u_end3,v_end3,end_p3,num_u3,num_v3,x_g3,y_g3,z_g3,var_cp_g3,h_g3,var_ws_g3,k_g3_u,k_g3_v,var_ks_g3_u,var_ks_g3_v,surf3);%construct the deformed patch1
end
if isempty(x3)
    srfo_4 = surf4;%if there is no variables for deformation surface, use the original surface
else
    [x_g4,y_g4,z_g4,h_g4,k_g4_u,k_g4_v,x_p4,y_p4,h_p4,k_p4]=paraSplit(x3,var_cp_g4,var_ws_g4,var_ks_g4_u,var_ks_g4_v,var_cp_p4,var_ws_p4,var_ks_p4);
    srfo_4= deformSurf(start_p4,u_end4,v_end4,end_p4,num_u4,num_v4,x_g4,y_g4,z_g4,var_cp_g4,h_g4,var_ws_g4,k_g4_u,k_g4_v,var_ks_g4_u,var_ks_g4_v,surf4);% construct the deformed patch2
end
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
u1_min = 0.0; u1_max =1.0;u2_min = 0.0; u2_max = 1.0;
% construct the points of the deformed patches
sum1=zeros;sum2=zeros;sum3=zeros;sum4=zeros;
if obj.num_trim1 ~= 0
    for i=1:obj.num_trim1
        pp = u1_min+(i-1)*(u1_max-u1_min)/(obj.num_trim1-1);
        %pt1 = trimCeval(tLoop1{num_Lp1_cstr},temp);
        obj.after1_trim1(i,:) =  ptOntrimCrv(srfo_1,trimCrv_p1.ini_value,pp);
        temp1 = (obj.after1_trim1(i,:)-obj.before1_trim1(i,:)).^2;% construct the points on trimming curve1 for optimization
        sum1 = sum1 +temp1;
    end
end
if obj.num_trim2 ~= 0
    for i=1:obj.num_trim2
        pp = u2_min+(i-1)*(u2_max-u2_min)/(obj.num_trim2-1);
       % pt2 = trimCeval(tLoop2{num_Lp2_cstr},pp);
        obj.after1_trim2(i,:) =  ptOntrimCrv(srfo_2,trimCrv_p2.ini_value,pp);
        temp2 = (obj.after1_trim2(i,:)-obj.before1_trim2(i,:)).^2;% construct the points on trimming curve2 for optimization
        sum2 = sum2+temp2;
    end
end
if obj.num_patch1 ~= 0
    for i=1:obj.num_patch1
        pp = obj.bp1{i};
        obj.after_patch1(i,:) =  nrbeval(srfo_1,{pp(1),pp(2)});
        temp1 = (obj.after_patch1(i,:)-obj.before_patch1(i,:)).^2;% construct the points on trimming curve1 for optimization
        sum1 = sum1 +temp1;
    end
end
if obj.num_patch2 ~= 0
    for i=1:obj.num_patch2
        pp = obj.bp2{i};
        obj.after_patch2(i,:) =  nrbeval(srfo_2,{pp(1),pp(2)});
        temp2 = (obj.after_patch2(i,:)-obj.before_patch2(i,:)).^2;% construct the points on trimming curve2 for optimization
        sum2 = sum2+temp2;
    end
end
if obj.num_patch3 ~= 0
    for i=1:obj.num_patch3
        pp = obj.bp3{i};
        obj.after_patch3(i,:) =  nrbeval(srfo_3,{pp(1),pp(2)});
        temp1 = (obj.after_patch3(i,:)-obj.before_patch3(i,:)).^2;% construct the points on trimming curve1 for optimization
        sum3 = sum3 +temp1;
    end
end
if obj.num_patch4 ~= 0
    for i=1:obj.num_patch4
        pp = obj.bp4{i};
        obj.after_patch4(i,:) =  nrbeval(srfo_4,{pp(1),pp(2)});
        temp2 = (obj.after_patch4(i,:)-obj.before_patch4(i,:)).^2;% construct the points on trimming curve2 for optimization
        sum4 = sum4+temp2;
    end
end

%% sumerize the total number
 sum=sum1+sum2+sum3+sum4;
 f= sum(1)+sum(2)+sum(3);
disp(sprintf('the objective function is %0.5f',f));
end

