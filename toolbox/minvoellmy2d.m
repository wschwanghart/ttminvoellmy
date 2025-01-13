classdef minvoellmy2d < handle
    %MINVOELLMY2D Version 2, 2024-02-12, Copyright (C) 2023-24 Stefan Hergarten 
    %This program is free software; you can redistribute it and/or modify it
    %under the terms of the GNU General Public License as published by the
    %Free Software Foundation; either version 3 of the License, or (at your
    %option) any later version.
    %This program is distributed in the hope that it will be useful, but
    %WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    %or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    %for more details. To view a copy of this license, visit
    %https://www.gnu.org/licenses/gpl-3.0.txt or write to the Free Software
    %Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 
    %For any questions concerning the codes and bug reports, please contact the
    %developer at stefan.hergarten@geologie.uni-freiburg.de. 
      
    properties
        b %Bed elevation of the nodes 
        h %Thickness of the mobile layer 
        uh %Momentum in x-direction
        vh %Momentum in y-direction
        wh %Momentum in z-direction
        dx %Grid spacing in x-direction
        dy %Grid spacing in y-direction
        dbdx %Gradient of the bed in x-direction
        dbdy %Gradient of the bed in y-direction
        cosbeta %Cosine of the slope angle
        stat %Flow status of each node (0 = not moving, 1 = Coulomb friction, 2 = Voellmy friction)
        mu %Coulomb coefficient of friction
        xi %Voellmy friction parameter
        vc %Crossover velocity at h = 1; vc <= 0 switches to the conventional Voellmy rheology
        g %Gravity (default = 9.81)
        hmin %Minimum thickness above which the layer can move (default = 0)
        dmin %Lower limit for the denominator if the original expression for the pressure shall be used (default = 0 = modified pressure)
        cent %Take into account the effect of centripetal acceleration on friction (default = true)
    end
    
    methods
        function o = minvoellmy2d(b,h,dx,dy,mu,xi,vc,hmin,dmin,cent,g)
            %MINVOELLMY2D Construct an instance of this class 
            %   All parameters are the same as the respective members of the class.
            %   b, h: 2D arrays of the same size.
            %   mu, xi, vc: scalar values or 2D arrays
            %   hmin, dmin: positive scalar values (default = 0)
            %   cent: logical (default = true)
            %   g: scalar value (default = 9.81)
                        
            o.b = b;
            o.h = h;
            o.dx = dx;
            o.dy = dy;
            % Central difference quotient
            o.dbdx = (o.b(:,[2:end,end-1])-o.b(:,[2,1:end-1]))/(2*o.dx);
            o.dbdy = (o.b([2:end,end-1],:)-o.b([2,1:end-1],:))/(2*o.dy);
            o.cosbeta = 1./sqrt(1+o.dbdx.^2+o.dbdy.^2);
            o.uh = zeros(size(o.h));
            o.vh = zeros(size(o.h));
            o.wh = zeros(size(o.h));
            o.mu = mu;
            o.xi = xi;
            o.vc = vc;
            if exist('hmin','var')
                o.hmin = hmin;
            else
                o.hmin = 0;
            end
            if exist('dmin','var')
                o.dmin = dmin;
            else
                o.dmin = 0;
            end
            if exist('cent','var')
                o.cent = cent;
            else
                o.cent = true;
            end
            if exist('g','var')
                o.g = g;
            else
                o.g = 9.81;
            end
        end
        
        function dt = step(o,dt,cfl)
            %STEP Forward time step
            %   The forward time step updates the members h, uh, vh, wh,
            %   and stat. If cfl is set, dt is automatically reduced to cfl
            %   times the limit defined by the CFL criterion. The used
            %   value of dt is returned.
            
            % Find rectangle around active nodes
            mask = o.h > o.hmin;
            cols = find(any(mask));
            cols = max(cols(1)-2,1):min(cols(end)+2,size(mask,2));
            rows = find(any(mask,2));
            rows = max(rows(1)-2,1):min(rows(end)+2,size(mask,1));
                        
            % Define data restricted to rectangle
            if numel(o.mu) > 1
                mu = o.mu(rows,cols);
            else
                mu = o.mu;
            end
            if numel(o.xi) > 1
                xi = o.xi(rows,cols);
            else
                xi = o.xi;
            end
            if numel(o.vc) > 1
                vc = o.vc(rows,cols);
            else
                vc = o.vc;
            end
            dbdx = o.dbdx(rows,cols);
            dbdy = o.dbdy(rows,cols);
            cosbeta = o.cosbeta(rows,cols);
            h = o.h(rows,cols);
            uh = o.uh(rows,cols);
            vh = o.vh(rows,cols);
            wh = o.wh(rows,cols);
            
            % x-velocity in the middle between the nodes
            u = uh./max(h,1e-10);
            u = (u(:,1:end-1)+u(:,2:end))/2;
            u(:,end+1) = 0;
            % y-velocity in the middle between the nodes
            v = vh./max(h,1e-10);
            v = (v(1:end-1,:)+v(2:end,:))/2;
            v(end+1,:) = 0;
                        
            % CFL criterion
            if nargin == 3; % exist('cfl','var')
                dtinv = max(abs(u)/o.dx+abs(v)/o.dy,[],'all')/cfl;
                if dtinv>1/dt
                    dt = 1/dtinv;
                end
            end
            
            % Upstream index for point in the middle between the nodes in x-direction 
            [ju,iu] = meshgrid(1:size(h,2),1:size(h,1));
            ju(u<0) = ju(u<0)+1;
            iu = sub2ind(size(h),iu,ju);
            % Upstream index for point in the middle between the nodes in y-direction 
            [jv,iv] = meshgrid(1:size(h,2),1:size(h,1));
            iv(v<0) = iv(v<0)+1;
            iv = sub2ind(size(h),iv,jv);
                        
            % Advection of h
            h = h-dt*(diff([zeros(size(h,1),1),h(iu).*u],1,2)/o.dx ...
                     +diff([zeros(1,size(h,2));h(iv).*v])/o.dy);
            % Advection of uh
            uh = uh-dt*(diff([zeros(size(h,1),1),uh(iu).*u],1,2)/o.dx ...
                       +diff([zeros(1,size(h,2));uh(iv).*v])/o.dy);
            % Advection of vh
            vh = vh-dt*(diff([zeros(size(h,1),1),vh(iu).*u],1,2)/o.dx ...
                       +diff([zeros(1,size(h,2));vh(iv).*v])/o.dy); 
            % Advection of wh
            wh = wh-dt*(diff([zeros(size(h,1),1),wh(iu).*u],1,2)/o.dx ...
                       +diff([zeros(1,size(h,2));wh(iv).*v])/o.dy);          
            % Clean up
            clear u v ju iu jv iv
            
            % Momentum
            mo = sqrt(uh.^2+vh.^2+wh.^2);
            % Product of centripetal acceleration, thickness, and dt
            hcdt = (uh.*dbdx+vh.*dbdy-wh).*cosbeta;
            % Make momentum parallel to bed
            uh = uh-hcdt.*dbdx.*cosbeta;
            vh = vh-hcdt.*dbdy.*cosbeta;
            wh = wh+hcdt.*cosbeta;
            % Rescale momentum to original value
            monew = sqrt(uh.^2+vh.^2+wh.^2);
            mask = monew>0;
            f = mo(mask)./monew(mask);
            uh(mask) = uh(mask).*f;
            vh(mask) = vh(mask).*f;
            wh(mask) = wh(mask).*f;
            
            % Gradient of surface height in x-direction
            dsdx = [ zeros(size(h,1),1), diff(o.b(rows,cols)+h,1,2)/o.dx, zeros(size(h,1),1) ];
            % Average with the right-hand and left-hand neighbor for
            % thickness-weighted central difference quotient
            hr = [ (h(:,1:end-1)+h(:,2:end))/2, zeros(size(h,1),1) ] + 1e-10;
            hl = circshift(hr,[0 1]);
            % Find local maxima in surface
            locmax = dsdx(:,1:end-1) > 0 & dsdx(:,2:end) < 0;
            % Suppress either hl or hr for local maxima in order to
            % switch to one-sided difference quotient
            mask = dsdx(:,1:end-1) < -dsdx(:,2:end);
            hl(locmax&mask) = 1e-10;
            hr(locmax&~mask) = 1e-10;
            dsdx = (dsdx(:,1:end-1).*hl+dsdx(:,2:end).*hr)./(hl+hr);
            
            % Gradient of surface height in y-direction
            dsdy = [ zeros(1,size(h,2)); diff(o.b(rows,cols)+h)/o.dy; zeros(1,size(h,2)) ];
            % Average with the right-hand and left-hand neighbor for
            % thickness-weighted central difference quotient
            hr = [ (h(1:end-1,:)+h(2:end,:))/2; zeros(1,size(h,2)) ] + 1e-10;
            hl = circshift(hr,[1 0]);
            % Find local maxima in surface
            locmax = dsdy(1:end-1,:) > 0 & dsdy(2:end,:) < 0;
            % Suppress either hl or hr for local maxima in order to
            % switch to one-sided difference quotient
            mask = dsdy(1:end-1,:) < -dsdy(2:end,:);
            hl(locmax&mask) = 1e-10;
            hr(locmax&~mask) = 1e-10;
            dsdy = (dsdy(1:end-1,:).*hl+dsdy(2:end,:).*hr)./(hl+hr);
            % Clean up           
            clear hr hl locmax mask   
             
            % Dot product of the gradients
            s = dsdx.*dbdx+dsdy.*dbdy;
            if o.dmin > 0
                % Original expression for pressure (not useful in general)
                p = o.g*h./max(1+s,o.dmin);
            else
                % Recommended expression for pressure
                p = o.g*h.*cosbeta.^2;
            end 
                                
            % Accelerate without friction
            uh = uh-dt*p.*dsdx;
            vh = vh-dt*p.*dsdy;
            wh = wh-dt*p.*s;
            % Clean up
            clear dsdx dsdy s
            
            % Find points with Voellmy friction
            mask = h > o.hmin & mo >= vc.*h.*(h.*cosbeta).^(1/3);
            
            mo = sqrt(uh.^2+vh.^2+wh.^2);
            monew = mo;
            % Apply Voellmy friction 
            f = xi.*h.^2.*cosbeta./(2*o.g*dt);
            f = f(mask);
            monew(mask) = max(sqrt(f.^2+2*f.*mo(mask))-f,0);
            
            if o.cent
                % Increase p by centrifugal acceleration
                p = max(p+hcdt.*cosbeta/dt,0);
            end 
            stat = 1+((vc>0)&mask)+(vc<=0)*(monew.^2>mu.*p.*xi.*h.^2./o.g);
                        
            % Find points with Coulomb friction
            mask = ~mask | vc <= 0;
               
            % Maximum momentum that can be consumed by friction
            maxf = mu.*p./cosbeta*dt;
            % Apply maximum Coulomb friction
            monew(mask) = monew(mask)-maxf(mask);
            
            % Find points that get stuck
            mask = h <= o.hmin | monew <= 0;
            monew(mask) = 0;
            stat(mask) = 0;
               
            % Rescale momentum and paste into orginal data
            f = monew./max(mo,1e-10);
            o.uh(rows,cols) = uh.*f;
            o.vh(rows,cols) = vh.*f;
            o.wh(rows,cols) = wh.*f;
            
            % Paste rest into original data
            o.h(rows,cols) = h;
            o.stat = zeros(size(o.b));
            o.stat(rows,cols) = stat;
        end  
    end
end