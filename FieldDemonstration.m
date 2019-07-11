#script to show a sample vector field created by two solenoids

 clf;
 colormap ("default");
 A = dlmread("Data.csv");

x = A(:,1)
y = A(:,2)
z = A(:,3)
 
u = A(:,4)
v = A(:,5)
w = A(:,6)

quiver3(x,y,z,u,v,w)

axis([-0.1 0.15 0.1 0.3 .04 .16])


 