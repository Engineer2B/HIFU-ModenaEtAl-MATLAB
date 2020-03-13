function son  = Reflection_Refraction(ray, object)
%Reflection_Reflaction Calculates the rays generated from a ray
%  Different Cases will be trated
son(15).xi=0;
%mat=Define_material();
cnt=1;

I_limit=0.00001;
ray.xi=1;
control=ray;
fr=1;

while isempty(control(1).xi)~=1 %if the control vector contains sons which can generate other sons
    sn=1; %contator for the control vector
    processing=control; %process the rays in the vector sons -> save them in processing
    control=[]; %clean control
    control(20).xi=[]; %initialize control
    [~,sizeP]=size(processing);
    fr=fr+1;
    
    for kk=1:sizeP %let's check them!
        %find the first zero of the path
        idx = find(processing(kk).path(1,:)==0, 1, 'first');
        
        switch (processing(kk).next_kind)
            
            case 0 % the next material is a soft tissue
                
                if object(processing(kk).actual_object).kind == 0 % ray is a soft tissue  %SOFT -> SOFT
                    %nn= check_sign(ray.Vray/norm(ray.Vray), object(ray.next_object).interface.normal); %check the normal of the plane
                    nn=processing(kk).nn; %the normal in the end point of the ray
                    [refr, v_out]=refraction0(processing(kk).Vray,nn,object(processing(kk).actual_object).c* 1e-2,object(processing(kk).next_object).c* 1e-2);%refractiononly long rays
                    
                    if processing(kk).actual_object==1 %ray is in the lossless
                        T=1; % everything is transmitted
                        R=0;
                    else
                        [T,R]=CoefficientLL(processing(kk).Vray,nn,object(processing(kk).actual_object).c* 1e-2,object(processing(kk).next_object).c* 1e-2,...%coefficient between two soft
                            object(processing(kk).actual_object).density* 1e6,object(processing(kk).next_object).density* 1e6);
                        
                    end
                    
                    if processing(kk).IF*T> I_limit % check if the refracted rays has enough initial energy
                        
                        if refr==1 && processing(kk).next_object ~= 1%refraction possible and the next mat is not lossless
                            son(cnt).actual_kind=0;
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0= processing(kk).IF*T;
                            son(cnt).Vray=v_out/norm(v_out);
                            if son(cnt).Vray(1)>0
                                son(cnt).direction=1; %going towards the focal point
                            else
                                son(cnt).direction=-1;
                            end
                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);  %calculate the intersection with all the objects
                            son(cnt).end= inters_point.pt; %point of the intersection
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            son(cnt).actual_object=processing(kk).next_object; %save the object in which the ray is
                            
                            
                            son(cnt).next_object=inters_point.index_obj; %intersection gives the object which the interface is belonging
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            
                            s2='final';
                            if strcmp(object(inters_point.index_obj).name_material,s2)%if the rays are at the end of the configuration
                                son(cnt).IF=0; % Set to zero so I don't have to treat them again
                            else
                                son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuation*inters_point.d);
                            end
                            son(cnt).phase_initial= processing(kk).phase_final; %the phase initial is the same of the father
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).k * inters_point.d ; %phase final
                            %son(cnt).previous_object=ray.actual_object; %really need it?
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object; %add the path to the ray
                            son(cnt).path=path_r;
                            if son(cnt).IF> I_limit %if the Intensity final of the sons can be able to generate other sons
                                %       if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in control
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn; % to be sure I am not quitting the while loop
                                sn=sn+1; %increase the sn number
                                %   end
                            end
                            if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                                ooooooooooo=0
                            end
                            cnt=cnt+1;
                            
                        end
                      
                    end
                    
                    if (R~=0 && processing(kk).actual_object ~= 1) %if R is not zero -> I have reflection!
                        
                        if processing(kk).IF*R> I_limit % check if the reflected ray has enough I0
                            v_out=reflection(processing(kk).Vray,nn);
                            son(cnt).Vray=v_out/norm(v_out);
                            if son(cnt).Vray(1)>0
                                son(cnt).direction=1; %going towards the focal point
                            else
                                son(cnt).direction=-1;
                            end
                            son(cnt).actual_kind=0;
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0= processing(kk).IF * R;

                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                            son(cnt).end= inters_point.pt;
                            son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                            son(cnt).next_object=inters_point.index_obj;
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            son(cnt).phase_initial= processing(kk).phase_final;
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).k * inters_point.d ;
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuation*inters_point.d);
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object;
                            son(cnt).path=path_r;
                            if son(cnt).IF> I_limit
                                % if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in the big Raylist
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn;
                                sn=sn+1;
                                % end
                            end
                            cnt=cnt+1;
                            
                        end
                        if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                            ooooooooooo=0
                        end
                    end
                    
                else % ray is in bone BONE -> MARROW
                    
                    
                    nn=processing(kk).nn;
                    alpha_in=acos(abs(nn'*processing(kk).Vray));
                    rho_liquid    = object(processing(kk).next_object).density*1e6;
                    rho_solid     = object(processing(kk).actual_object).density*1e6;
                    c_long_liquid = object(processing(kk).next_object).c*1e-2;
                    c_long_solid  = object(processing(kk).actual_object).clong*1e-2;
                    c_shear_solid = object(processing(kk).actual_object).cshear*1e-2;
                    
                    %longitudinal wave -> reflected shear, reflected long, refracted long
                    if processing(kk).shear==0
                        
                        [Refl_long,Refl_shear,Transm_long]= B2MReflTransmCoef_LongInp(alpha_in,rho_liquid,rho_solid,c_long_liquid,...
                            c_long_solid,c_shear_solid);
                        
                        intens_long_reflec= processing(kk).IF*Refl_long;
                        
                        %REFLECTED LONGITUDINAL RAY
                        
                        if  intens_long_reflec>I_limit
                            v_out= reflection(processing(kk).Vray,nn);
                            son(cnt).Vray=v_out/norm(v_out);
                            if son(cnt).Vray(1)>0
                                son(cnt).direction=1; %going towards the focal point
                            else
                                son(cnt).direction=-1;
                            end
                            son(cnt).actual_kind=1;
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0=  intens_long_reflec;
                            
                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                            son(cnt).end= inters_point.pt;
                            son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                            son(cnt).next_object=inters_point.index_obj;
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            son(cnt).phase_initial= processing(kk).phase_final;
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).kl * inters_point.d ;
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuationl*inters_point.d);
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            son(cnt).shear=0;
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object;
                            son(cnt).path=path_r;
                            if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                                son(cnt).next_object=3;
                                son(cnt).next_kind=0;
                            end
                            if son(cnt).IF> I_limit
                                % if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in the big Raylist
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn;
                                sn=sn+1;
                                % end
                                
                            end
                            cnt=cnt+1;
                        end
                        
                        %   REFLECTED SHEAR RAY
                        intens_shear_reflec= processing(kk).IF*Refl_shear;
                        
                        if intens_shear_reflec>I_limit
                            [refl,v_out,poldir]= reflection2(processing(kk).Vray,nn,c_long_solid,c_shear_solid);
                            
                            son(cnt).Vray=v_out/norm(v_out);
                            if son(cnt).Vray(1)>0
                                son(cnt).direction=1; %going towards the focal point
                            else
                                son(cnt).direction=-1;
                            end
                            son(cnt).actual_kind=1;
                            son(cnt).polarization=poldir; %it's a shear save the polarization direction
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0= intens_shear_reflec;
                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                            son(cnt).end= inters_point.pt;
                            son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                            son(cnt).next_object=inters_point.index_obj;
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            son(cnt).phase_initial= processing(kk).phase_final;
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).ks * inters_point.d ;
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuations*inters_point.d);
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            son(cnt).shear=1; %It's a shear
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object;
                            son(cnt).path=path_r;
                            if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                                 son(cnt).next_object=3;
                                 son(cnt).next_kind=0;
                            end
                            
                            
                            if son(cnt).IF> I_limit
                                % if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in the big Raylist
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn;
                                sn=sn+1;
                                % end
                            end
                            cnt=cnt+1;
                        end
                        
                        
                        %REFRACTED LONGITUDINAL RAY
                        intens_long_refract=processing(kk).IF*Transm_long;
                        if intens_long_refract>I_limit
                            [refr, v_out]=refraction(processing(kk).Vray,nn,c_long_solid,c_long_liquid);
                            if refr==1
                                son(cnt).actual_kind=0;
                                son(cnt).start=processing(kk).end;
                                son(cnt).I0= intens_long_refract;
                                son(cnt).Vray=v_out/norm(v_out);
                                if son(cnt).Vray(1)>0
                                    son(cnt).direction=1; %going towards the focal point
                                else
                                    son(cnt).direction=-1;
                                end
                                inters_point= Intersection(son(cnt),object,processing(kk).next_object);  %calculate the intersection with all the objects
                                son(cnt).end= inters_point.pt; %point of the intersection
                                son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                                son(cnt).actual_object=processing(kk).next_object; %save the object in which the ray is                 
                                son(cnt).next_object=inters_point.index_obj; %intersection gives the object which the interface is belonging
                                son(cnt).next_kind=object(son(cnt).next_object).kind;
                                
                                s2='final';
                                if strcmp(object(inters_point.index_obj).name_material,s2)%if the rays are at the end of the configuration
                                    son(cnt).IF=0; % Set to zero so I don't have to treat them again
                                else
                                    son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuation*inters_point.d);
                                end
                                son(cnt).phase_initial= processing(kk).phase_final; %the phase initial is the same of the father
                                son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).k * inters_point.d ; %phase final
                                %son(cnt).previous_object=ray.actual_object; %really need it?
                                path_r= processing(kk).path;
                                path_r(1, idx)=son(cnt).actual_object; %add the path to the ray
                                son(cnt).path=path_r;
                                
                                if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                                ooooooooooo=0
                            end
                                
                                
                                
                                if son(cnt).IF> I_limit %if the Intensity final of the sons can be able to generate other sons
                                    %       if son(cnt).actual_object~=1
                                    fn=fieldnames(son);% fields of son
                                    [s,~]=size(fn);
                                    for jj=1:s  %copy the rays in control
                                        [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                    end
                                    control(sn).xi=sn; % to be sure I am not quitting the while loop
                                    sn=sn+1; %increase the sn number
                                    %   end
                                end
                                cnt=cnt+1;
                            end
                            
                        end
                        
                        
                        
                    else % it's a shear wave. Shear: reflected shear horiz pol, reflected shear vert pol, reflected long, refracted long
                        nn=processing(kk).nn;
                        [v_out, verticalPolarisationIn,verticalPolarisationOut,horizontalPolarisation]=reflection3(processing(kk).Vray,nn);
                        %alpha_in=acos(abs(nn'*processing(kk).Vray));
                        pol_dir_in=processing(kk).polarization;
                        cosvert=pol_dir_in'*verticalPolarisationIn;
                        coshor=pol_dir_in'*horizontalPolarisation;
                        intensity_horcomp= processing(kk).IF* coshor^2;
                        if intensity_horcomp > I_limit % SHEAR HORIZONTAL REFLECTED
                            
                            son(cnt).Vray=v_out/norm(v_out);
                            if son(cnt).Vray(1)>0
                                son(cnt).direction=1; %going towards the focal point
                            else
                                son(cnt).direction=-1;
                            end
                            son(cnt).actual_kind=1;
                            son(cnt).polarization=horizontalPolarisation; %it's a shear save the polarization direction
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0= intensity_horcomp;
                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                            son(cnt).end= inters_point.pt;
                            son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                            son(cnt).next_object=inters_point.index_obj;
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            son(cnt).phase_initial= processing(kk).phase_final;
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).ks * inters_point.d ;
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuations*inters_point.d);
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            son(cnt).shear=1; %It's a shear
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object;
                            son(cnt).path=path_r;
                            if son(cnt).IF> I_limit
                                % if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in the big Raylist
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn;
                                sn=sn+1;
                                % end
                                
                                
                            end
                            cnt=cnt+1;
                        end
                       
                        alpha_in=acos(abs(nn'*processing(kk).Vray));
                        rho_liquid    = object(processing(kk).next_object).density*1e6;
                        rho_solid     = object(processing(kk).actual_object).density*1e6;
                        c_long_liquid = object(processing(kk).next_object).c*1e-2;
                        c_long_solid  = object(processing(kk).actual_object).clong*1e-2;
                        c_shear_solid = object(processing(kk).actual_object).cshear*1e-2;
                        [Refl_long,Refl_shear,Transm_long]= B2MReflTransmCoef_ShearInp(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid);
                        intensity_vertcomp=processing(kk).IF* cosvert^2;
                        intensity_vert_shear_reflected=Refl_shear* intensity_vertcomp;
                        
                        if intensity_vert_shear_reflected> I_limit
                            son(cnt).Vray=v_out/norm(v_out); %same V_out reflection3
                            if son(cnt).Vray(1)>0
                                son(cnt).direction=1; %going towards the focal point
                            else
                                son(cnt).direction=-1;
                            end
                            son(cnt).actual_kind=1;
                            son(cnt).polarization=verticalPolarisationOut; %it's a shear save the polarization direction
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0= intensity_vert_shear_reflected;
                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                            son(cnt).end= inters_point.pt;
                            son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                            son(cnt).next_object=inters_point.index_obj;
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            son(cnt).phase_initial= processing(kk).phase_final;
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).ks * inters_point.d ;
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuations*inters_point.d);
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            son(cnt).shear=1; %It's a shear
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object;
                            son(cnt).path=path_r;
                             if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                            ooooooooooo
                        end
                            if son(cnt).IF> I_limit
                                % if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in the big Raylist
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn;
                                sn=sn+1;
                                % end
                                
                                
                            end
                            cnt=cnt+1;
                        end
                        % LONGITUDINAL REFLECTED
                        
                        [longrefl_possible,v_out, poldir]=reflection2(processing(kk).Vray,nn,c_shear_solid,c_long_solid);
                         intensity_longit_reflected=Refl_long* intensity_vertcomp; 
                         if longrefl_possible==1 && intensity_longit_reflected>I_limit
                             son(cnt).Vray=v_out/norm(v_out);
                             if son(cnt).Vray(1)>0
                                 son(cnt).direction=1; %going towards the focal point
                             else
                                 son(cnt).direction=-1;
                             end
                             son(cnt).actual_kind=1;
                            son(cnt).start=processing(kk).end;
                            son(cnt).I0=  intensity_longit_reflected;
                            
                            inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                            son(cnt).end= inters_point.pt;
                            son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                            son(cnt).next_object=inters_point.index_obj;
                            son(cnt).next_kind=object(son(cnt).next_object).kind;
                            son(cnt).phase_initial= processing(kk).phase_final;
                            son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).kl * inters_point.d ;
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuationl*inters_point.d);
                            son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                            son(cnt).shear=0;
                            path_r= processing(kk).path;
                            path_r(1, idx)=son(cnt).actual_object;
                            son(cnt).path=path_r;
                             if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                            ooooooooooo
                        end
                            
                            
                            if son(cnt).IF> I_limit
                                % if son(cnt).actual_object~=1
                                fn=fieldnames(son);% fields of son
                                [s,~]=size(fn);
                                for jj=1:s  %copy the rays in the big Raylist
                                    [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                end
                                control(sn).xi=sn;
                                sn=sn+1;
                                % end
                                
                            end
                            cnt=cnt+1;
                         end
                         
                         %LONGITUDINAL REFRACTED
                         [refr, v_out,poldir]=refraction(processing(kk).Vray,nn,c_shear_solid,c_long_liquid);
                         intensity_longit_refracted=Transm_long* intensity_vertcomp;
                         
                         if intensity_longit_refracted>I_limit
                             son(cnt).actual_kind=0;
                             son(cnt).start=processing(kk).end;
                             son(cnt).I0= intensity_longit_refracted;
                             son(cnt).Vray=v_out/norm(v_out);
                             if son(cnt).Vray(1)>0
                                 son(cnt).direction=1; %going towards the focal point
                             else
                                 son(cnt).direction=-1;
                             end
                             inters_point= Intersection(son(cnt),object,processing(kk).next_object);  %calculate the intersection with all the objects
                             son(cnt).end= inters_point.pt; %point of the intersection
                             son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                             son(cnt).actual_object=processing(kk).next_object; %save the object in which the ray is
                             son(cnt).next_object=inters_point.index_obj; %intersection gives the object which the interface is belonging
                             son(cnt).next_kind=object(son(cnt).next_object).kind;
                             
                             s2='final';
                             if strcmp(object(inters_point.index_obj).name_material,s2)%if the rays are at the end of the configuration
                                 son(cnt).IF=0; % Set to zero so I don't have to treat them again
                             else
                                 son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuation*inters_point.d);
                             end
                             son(cnt).phase_initial= processing(kk).phase_final; %the phase initial is the same of the father
                             son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).k * inters_point.d ; %phase final
                             %son(cnt).previous_object=ray.actual_object; %really need it?
                             path_r= processing(kk).path;
                             path_r(1, idx)=son(cnt).actual_object; %add the path to the ray
                             son(cnt).path=path_r;
                             
                             
                              if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                            ooooooooooo
                        end
                             
                             if son(cnt).IF> I_limit %if the Intensity final of the sons can be able to generate other sons
                                 %       if son(cnt).actual_object~=1
                                 fn=fieldnames(son);% fields of son
                                 [s,~]=size(fn);
                                 for jj=1:s  %copy the rays in control
                                     [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                                 end
                                 control(sn).xi=sn; % to be sure I am not quitting the while loop
                                 sn=sn+1; %increase the sn number
                                 %   end
                             end
                             cnt=cnt+1;
                             
                             
                         end
                       
                         
                    end %if processing(kk).shear==0
                    
                    
                end
                
                
                
            case 1 % next material is bone. SOFT->BONE
                
                nn=processing(kk).nn;
                alpha_in=acos(abs(nn'*processing(kk).Vray));
                rho_liquid    = object(processing(kk).actual_object).density*1e6;
                rho_solid     = object(processing(kk).next_object).density*1e6;
                c_long_liquid = object(processing(kk).actual_object).c*1e-2;
                c_long_solid  = object(processing(kk).next_object).clong*1e-2;
                c_shear_solid = object(processing(kk).next_object).cshear*1e-2;
                %calculation coefficient
                if isempty (c_long_liquid)
                    dhdhdhdhdhdhdh
                end
                
                [Refl, Transm_long,Transm_shear]= M2BReflTransmCoef(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid);
                
                % REFLECTED longitudinal ray in soft tissue
                
                intens_long_reflec= processing(kk).IF*Refl;% intensity of the reflected
                
                if intens_long_reflec>I_limit %if the reflected has enough initial energy
                    
                    v_out=reflection(processing(kk).Vray,nn);
                    son(cnt).actual_kind=0;
                    son(cnt).start=processing(kk).end;
                    son(cnt).I0= intens_long_reflec;
                    son(cnt).Vray=v_out/norm(v_out);
                    if son(cnt).Vray(1)>0
                        son(cnt).direction=1; %going towards the focal point
                    else
                        son(cnt).direction=-1;
                    end
                    inters_point= Intersection(son(cnt),object,processing(kk).next_object);
                    
                    son(cnt).end= inters_point.pt;
                    son(cnt).actual_object=processing(kk).actual_object; %the reflected is in the same material
                    son(cnt).next_object=inters_point.index_obj;
                    son(cnt).next_kind=object(son(cnt).next_object).kind;
                    son(cnt).phase_initial= processing(kk).phase_final;
                    son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).k * inters_point.d ;
                    son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuation*inters_point.d);
                    son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                    path_r= processing(kk).path;
                    
                    path_r(1, idx)=son(cnt).actual_object;
                    son(cnt).path=path_r;
                    
                    
                    if son(cnt).IF> I_limit
                        %                         % if son(cnt).actual_object~=1
                        fn=fieldnames(son);% fields of son
                        [s,~]=size(fn);
                        for jj=1:s  %copy the rays in the big Raylist
                            [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                        end
                        control(sn).xi=sn;
                        sn=sn+1;
                        %                         % end
                        
                    end
                    cnt=cnt+1;
                end %end of the reflected
                
                % REFRACTED longitudinal ray in BONE
                [longrefrac_possible,v_out]=refraction(processing(kk).Vray,nn, c_long_liquid,c_long_solid);
                intens_long_refrac= processing(kk).IF*Transm_long;
                if  longrefrac_possible==1 && intens_long_refrac>I_limit
                    
                    son(cnt).actual_kind=1;
                    son(cnt).start=processing(kk).end;
                    son(cnt).I0= intens_long_refrac;
                    son(cnt).Vray=v_out/norm(v_out);
                    if son(cnt).Vray(1)>0
                        son(cnt).direction=1; %going towards the focal point
                    else
                        son(cnt).direction=-1;
                    end
                    inters_point= Intersection(son(cnt),object,processing(kk).next_object);  %calculate the intersection with all the objects
                    
                    son(cnt).end= inters_point.pt; %point of the intersection
                    son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                    son(cnt).actual_object=processing(kk).next_object; %save the object in which the ray is
                    son(cnt).next_object=inters_point.index_obj; %intersection gives the object which the interface is belonging
                    
                    
                    son(cnt).next_kind=object(son(cnt).next_object).kind;
                    s2='final';
                    if strcmp(object(inters_point.index_obj).name_material,s2)%if the rays are at the end of the configuration
                        son(cnt).IF=0; % Set to zero so I don't have to treat them again
                    else
                        son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuationl*inters_point.d);
                    end
                    son(cnt).phase_initial= processing(kk).phase_final; %the phase initial is the same of the father
                    son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).kl * inters_point.d ; %phase final
                    %son(cnt).previous_object=ray.actual_object; %really need it?
                    path_r= processing(kk).path;
                    path_r(1, idx)=son(cnt).actual_object; %add the path to the ray
                    son(cnt).path=path_r;
                    son(cnt).shear=0; % I know it's a longitudinal
                    
                    if son(cnt).IF> I_limit %if the Intensity final of the sons can be able to generate other sons
                        %       if son(cnt).actual_object~=1
                        fn=fieldnames(son);% fields of son
                        [s,~]=size(fn);
                        for jj=1:s  %copy the rays in control
                            [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                        end
                        control(sn).xi=sn; % to be sure I am not quitting the while loop
                        sn=sn+1; %increase the sn number
                        %   end
                    end
                    cnt=cnt+1;
                    
                end
                % REFRACTED Shear ray in BONE
                
                [shearrefrac_possible,v_out,poldir]=refraction(processing(kk).Vray,nn,c_long_liquid,c_shear_solid);
                if shearrefrac_possible==1 % there is shear
                    
                    intens_shear_refrac=processing(kk).IF*Transm_shear;
                    
                    if intens_shear_refrac>I_limit %the intensity of the shear is enough
                        
                        son(cnt).actual_kind=1;
                        son(cnt).start=processing(kk).end;
                        son(cnt).I0= intens_shear_refrac;
                        son(cnt).Vray=v_out/norm(v_out);
                        son(cnt).polarization=poldir; %it's a shear save the polarization direction
                        if son(cnt).Vray(1)>0
                            son(cnt).direction=1; %going towards the focal point
                        else
                            son(cnt).direction=-1;
                        end
                        inters_point= Intersection(son(cnt),object,processing(kk).next_object);  %calculate the intersection with all the objects
                        son(cnt).end= inters_point.pt; %point of the intersection
                        son(cnt).nn=  inters_point.normal; % the normal at the end of the ray
                        son(cnt).actual_object=processing(kk).next_object; %save the object in which the ray is
                        son(cnt).next_object=inters_point.index_obj; %intersection gives the object which the interface is belonging
                        son(cnt).next_kind=object(son(cnt).next_object).kind;
                        s2='final';
                        if strcmp(object(inters_point.index_obj).name_material,s2)%if the rays are at the end of the configuration
                            son(cnt).IF=0; % Set to zero so I don't have to treat them again
                        else
                            son(cnt).IF=son(cnt).I0*exp(-2*object(son(cnt).actual_object).attenuations*inters_point.d);
                        end
                        son(cnt).phase_initial= processing(kk).phase_final; %the phase initial is the same of the father
                        son(cnt).phase_final= processing(kk).phase_final+ object(son(cnt).actual_object).ks * inters_point.d ; %phase final
                        %son(cnt).previous_object=ray.actual_object; %really need it?
                        path_r= processing(kk).path;
                        path_r(1, idx)=son(cnt).actual_object; %add the path to the ray
                        son(cnt).path=path_r;
                        son(cnt).shear=1; % I know it's a shear
                        if son(cnt).direction==1 & son(cnt).next_object== son(cnt).actual_object
                            ooooooooooo
                        end
                        
                       
                            
                        if son(cnt).IF> I_limit %if the Intensity final of the sons can be able to generate other sons
                            %       if son(cnt).actual_object~=1
                            fn=fieldnames(son);% fields of son
                            [s,~]=size(fn);
                            for jj=1:s  %copy the rays in control
                                [control(sn).(fn{jj})]=son(cnt).(fn{jj}); %copying the fields
                            end
                            control(sn).xi=sn; % to be sure I am not quitting the while loop
                            sn=sn+1; %increase the sn number
                            %   end
                        end
                        cnt=cnt+1;
                        
                        
                    end
                end
                
                
        end %switch
    end %sizeP
    
    if ~isempty(control(1).xi) %clean the control vector
        vec_mat={control.xi}';
        svid=cellfun(@isempty, vec_mat);
        control(svid)=[];
    end
    
    
    
    
end %end of while

% %clean the son vector
% if ~isempty(son(1).xi)
%     vec_mat={son.actual_object}';
%     svid=cellfun(@isempty, vec_mat);
%     son(svid)=[];
% end

end

