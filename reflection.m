function  v_out=reflection(v_in,n);
% reflection with incoming direction v_in, normal n , 
% v_in, n and v_out must be normalised to length
v_in=v_in/norm(v_in);
n=n/norm(n);
v_out=v_in-2*(v_in'*n)*n;



