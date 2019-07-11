/* Solenoid Analysis
 * Program simulates B field around two solenoids using biot savart law
 * Victor Guo
 */
#include<bits/stdc++.h>

using namespace std;
typedef long long ll;

//Create struct of vectors
struct vect{
public:
    double x,y,z; //3D vectors
    //Constructor
    vect(double x_in, double y_in, double z_in) {
        x = x_in; y = y_in; z = z_in;
    }
};
//Method to compute the cross product of two vectors
vect crossProduct(vect a, vect b){
    return vect(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y * b.x);
}
//Method to compute scalar product of vector and scalar
vect scalarProduct(vect a, double k){
    return vect(a.x*k, a.y*k, a.z*k);
}
//Method to add two vectors
vect addVector(vect a, vect b){
    return vect(a.x+b.x, a.y+b.y, a.z+b.z);
}
//Method to compute the magnitude of a vector
double magnitude(vect a){
    return pow(a.x*a.x + a.y*a.y + a.z*a.z, 0.5);
}
//Method to normalize a vector
vect normalize(vect a){
    return vect(a.x/magnitude(a), a.y/magnitude(a), a.z/magnitude(a));
}


//Variable Declarations
#define PI 3.14159265358
vector<pair<vect,vect> > magField; //array relating a position vector to a magnetic field strength
double h1 = 0.055, w1 = 0.055, l1 = 0.066; //dummy values
double I = 0.06, I2 = 0.06;
const double u = 4*PI*0.0000001; //permittivity constant
double h2 = 0.055, w2 = 0.055, l2 = 0.066; //dummy values

void input(){
    printf("This program computes the field between two rectangular cylindrical solenoids.\n");
    printf("Enter the dimensions for the first solenoid: \n");
    printf("Height: ");
    scanf("%lf", &h1);
    printf("Width: ");
    scanf("%lf", &w1);
    printf("Length: ");
    scanf("%lf", &l1);
    printf("Current: ");
    scanf("%lf", &I);

    printf("\n");

    printf("Enter the dimensions for the second solenoid: \n");
    printf("Height: ");
    scanf("%lf", &h2);
    printf("Width: ");
    scanf("%lf", &w2);
    printf("Length: ");
    scanf("%lf", &l2);
    printf("Current: ");
    scanf("%lf", &I2);
}

//Method to initialize experiment setup
void initialize(){
 for(int i = 0; i < 10; i++){
    for(int j = 0; j < 10; j++){
        for(int k = 0; k < 10; k++){
            magField.push_back(make_pair(vect(-0.15 + i*3.0/100.0, 3.0/200 + 3.0*j/100.0, 1.0/20 + 1.0*k/100.0), vect(0,0,0)));
        }
    }
 }
}

//Method to calculate magnetic field at position
vect calculateField(vect pos){
    vect B = vect(0,0,0);
    vect r = vect(0,0,0);
    vect L = vect(0,0,0);
    //square coil
    //xi, yi, zi instantaneous position
    double xi,yi,zi;
    //top
    zi = 0.1+h1;
    L = vect(w1/50, 0, 0);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 45; j++){
            xi = -w1/2.0 + 1.0/100*w1 + i/50.0*w1;
            yi = 0.1 + l1 * (j-1.0)/44.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //bottom
    zi = 0.1;
    L = vect(-w1/50, 0, 0);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 45; j++){
            xi = -w1/2.0 + 1.0/100*w1 + i/50.0*w1;
            yi = 0.1 + l1 * (j-1.0)/44.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //right
    xi = w1/2;
    L = vect(0,0,-h1/50);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 45; j++){
            zi = 0.1+(2*i+1.0)*h1/100.0;
            yi = 0.1 + l1 * (j-1.0)/44.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //left
    xi = -w1/2;
    L = vect(0,0,h1/50);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 45; j++){
            zi = 0.1+(2*i+1.0)*h1/100.0;
            yi = 0.1 + l1 * (j-1.0)/44.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //solenoid 2
    //top
    zi = 0.05+h2;
    L = vect(w2/50, 0, 0);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 90; j++){
            xi = -w2/2.0 + 1.0/100*w2 + i/50.0*w2;
            yi = 0.13 + l1+ l2 * (j-1.0)/89.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I2/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //bottom
    zi = 0.05;
    L = vect(-w2/50, 0, 0);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 90; j++){
            xi = -w2/2.0 + 1.0/100*w2 + i/50.0*w2;
            yi = 0.13 + l1+ l2 * (j-1.0)/89.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I2/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //right
    xi = w2/2;
    L = vect(0,0,-h2/50);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 90; j++){
            zi = 0.05+(2*i+1.0)*h1/100.0;
            yi = 0.13 + l1+ l2 * (j-1.0)/89.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I2/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }
    //left
    xi = -w2/2;
    L = vect(0,0,h2/50);
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 90; j++){
            zi = 0.05+(2*i+1.0)*h1/100.0;
            yi = 0.13 + l1+ l2 * (j-1.0)/89.0;
            r = addVector(pos, vect(-xi,-yi,-zi));
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I2/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }

    //Code for a secondary circular coil
    /*
    for(int i = 0; i < 50; i++){
        for(int j = 1; j <= 50; j++){
            yi = 0.1+l1+0.05+l2*(j-1.0)/50;
            xi = r2*cos(2*PI*i/50);
            zi = r2+h2+r2*sin(2*PI*i/50);
            r = addVector(pos,vect(-xi,-yi,-zi));
            L = vect(-sin(2*PI*i/50)*2*PI*r2/50.0,0,cos(2*PI*i/50)*2*PI*r2/50.0);
            B = addVector(B, scalarProduct(crossProduct(L,normalize(r)),u*I/(4.0 * PI * magnitude(r) * magnitude(r))));
        }
    }*/
    return B;
}
int main(){
    //input();
    initialize();
    for(int i = 0; i < magField.size(); i++){
        magField[i].second = calculateField(magField[i].first);
    }
    freopen ("Data.txt","w",stdout);
    for(auto val : magField){
        printf("%f %f %f %f %f %f\n", val.first.x, val.first.y, val.first.z, val.second.x*1000000.0, val.second.y*1000000.0, val.second.z*1000000.0); //output in microtesla
    }
    fclose (stdout);
	return 0;
}
