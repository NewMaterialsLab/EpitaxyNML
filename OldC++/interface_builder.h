 #include <vector>
 #include <iostream>
 #include <fstream>
 #include <string>
 #include <cmath>
 #include <cstdlib>
 #include <iomanip>
 #include <algorithm>
 using namespace std;

 const double maxMCIA = 100.0;
 const double maxAreaErr = 12.0;

 
 //Read Surface POSCAR
 std::vector<vector<double> > readSurf(string);
 //Calculate N1 and N2 for conmensurability
 std:: vector<vector<int> > calcAreasList(double, double);
 //Calculate Area for 2x2 matrix
 double calcArea(vector<vector<double> >);
 //Comapre superlattices N1,N2
 void compareSupperLattice(int, int);
 //Comapre superlattices N1,N2
 void readPlanes();
 //Print List
 void printList();
 //Sort List
 void sortList();
 //Create superslab matrix generator
 std::vector<vector<vector<int> > > matGen(int);

 // Unit surface cell for thin filmi and substrate POSCAR 
 vector<vector<double> > surfTF(2, vector<double>(2));
 vector<vector<double> > surfSub(2, vector<double>(2));

 // Planes input
 vector<vector<int> > planesTF;
 vector<vector<int> > planesSub;

 //Result list
 vector<vector<double> > listR;
