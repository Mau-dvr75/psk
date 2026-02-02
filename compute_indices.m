% calculates the indices at which objects (sources, receivers) reside,
% based on their location in true coordinates

function [coord_x_id, coord_z_id, nthings] = compute_indices(coord_x,coord_z,dh,nx,nz)

    nthings=length(coord_x);
    
    coord_x_id=zeros(1,nthings);
    coord_z_id=zeros(1,nthings);

    x=[0:nx-1]*dh;
    z=[-2:nz-1]*dh;

    for i=1:nthings
        coord_x_id(i)=min(find(min(abs(x-coord_x(i)))==abs(x-coord_x(i))));
        coord_z_id(i)=min(find(min(abs(z-coord_z(i)))==abs(z-coord_z(i))));
        %if i>1
         %   coord_x_id(i)=coord_x_id(i)+1;
        %end
        %[d,coord_x_id2(i)]=min(abs(x-coord_x(i)));
    end

end