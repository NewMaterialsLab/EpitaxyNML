 #include "interface_builder.h"
 
 vector<vector<double> > readSurf(string namefile){
   //cout << namefile<< std::endl;
   vector<vector<double> > kk(2, vector<double>(2));
   vector<vector<double> > surf3D(3, vector<double>(3));
   std::ifstream infile;
   infile.open(namefile.c_str());
   
   // Commented line
   string comments;
   getline(infile,comments);
   // Factor line
   double factor;
   infile >> factor; 
   // Lattice
   for(int i=0; i<3; i++){
     for(int j=0; j<3; j++){
       infile >>  surf3D[i][j];
     }
   }
   // Obtain surface lattice * factor
   for(int i=0; i<2; i++){
     for(int j=0; j<2; j++){
       kk[i][j] = surf3D[i][j]*factor;
     }
   }
   return kk; 
 } 

 vector<vector<int> > calcAreasList(double areaPrimTF, double areaPrimSub){
   double MCIA = 0.0;
   int it=1; 
   vector<vector<int> > listAreas;
   vector<int> dummyVec(2);
   while(MCIA < maxMCIA){
     dummyVec[0]=it;
     MCIA = areaPrimSub*it;
	 //cout << "it = " << it << "MCIA = " << MCIA << std::endl;
     int it2 = 1;
     double dummyArea = 0.0;
     while(dummyArea < maxMCIA){
       dummyVec[1]=it2;
       dummyArea = areaPrimTF*it2;
	 	//cout << "it2 = " << it2 << "dummyArea = " << dummyArea << std::endl;
       double difArea = 100*(dummyArea-MCIA)/dummyArea;
		//cout << "difArea = " << difArea <<std::endl; 
       if(abs(difArea)<maxAreaErr){
         //cout << "Matching Areas " << dummyVec[0] << " " << dummyVec[1] << std::endl;
		 //cout << "MATCH " << std::endl;
         listAreas.push_back(dummyVec);
       }
       it2 = it2 + 1;
     }
     it = it + 1;
   }   
   return listAreas; 
 }

 double calcArea(vector<vector<double> > lattice){
   double kk; 
   kk = abs((lattice[0][0]*lattice[1][1])-(lattice[1][0]*lattice[0][1]));
   return kk; 
 }

 void compareSupperLattice(int nSub, int nTF){
   vector<vector<vector<int> > > genMatSub, genMatTF;
   // Create superslabs generator matrixes
   genMatSub = matGen(nSub);
   cout << " Printing MatGen for Sub " << nSub << std::endl;
   for(int i=0; i<genMatSub.size(); i++){
	 cout << genMatSub[i][0][0] << " " << genMatSub[i][0][1] << " " << genMatSub[i][1][0] << " " << genMatSub[i][1][1] << std::endl; 
	}
   genMatTF =  matGen(nTF);
   cout << " Printing MatGen for TF " << nSub << std::endl;
   for(int i=0; i<genMatTF.size(); i++){
	 cout << genMatTF[i][0][0] << " " << genMatTF[i][0][1] << " " << genMatTF[i][1][0] << " " << genMatTF[i][1][1] << std::endl; 
	}

   //Create superslab for Sub
   vector<vector<double> >superSlabSub;
   for(int i=0; i<genMatSub.size(); i++){ 
     vector<double> kk(2);   
     vector<double> kk2(2);   
     //Generate new vectors
     kk[0]=genMatSub[i][0][0]*surfSub[0][0]+genMatSub[i][0][1]*surfSub[1][0];
     kk[1]=genMatSub[i][0][0]*surfSub[0][1]+genMatSub[i][0][1]*surfSub[1][1];
     kk2[0]=genMatSub[i][1][1]*surfSub[1][0];
     kk2[1]=genMatSub[i][1][1]*surfSub[1][1];
     // Compute modules and angles
     double modulkk, modulkk2, alpha;
     modulkk = sqrt(kk[0]*kk[0]+kk[1]*kk[1]);
     modulkk2 = sqrt(kk2[0]*kk2[0]+kk2[1]*kk2[1]);
     alpha = acos((kk[0]*kk2[0]+kk[1]*kk2[1])/(modulkk*modulkk2)); 
     vector<double> entry(3);
     if(modulkk<modulkk2){
       entry[0] = modulkk;
       entry[1] = modulkk2;
     }else{
       entry[0] = modulkk2;
       entry[1] = modulkk;
     }
     entry[2] = alpha;
     superSlabSub.push_back(entry);
   }

   //Create superslab for Sub
   vector<vector<double> >superSlabTF;
   for(int i=0; i<genMatTF.size(); i++){ 
     vector<double> kk(2);   
     vector<double> kk2(2);   
     //Generate new vectors
     kk[0]=genMatTF[i][0][0]*surfTF[0][0]+genMatTF[i][0][1]*surfTF[1][0];
     kk[1]=genMatTF[i][0][0]*surfTF[0][1]+genMatTF[i][0][1]*surfTF[1][1];
     kk2[0]=genMatTF[i][1][1]*surfTF[1][0];
     kk2[1]=genMatTF[i][1][1]*surfTF[1][1];
     // Compute modules and angles
     double modulkk, modulkk2, alpha;
     modulkk = sqrt(kk[0]*kk[0]+kk[1]*kk[1]);
     modulkk2 = sqrt(kk2[0]*kk2[0]+kk2[1]*kk2[1]);
     alpha = acos((kk[0]*kk2[0]+kk[1]*kk2[1])/(modulkk*modulkk2)); 
     vector<double> entry(3);
     if(modulkk<modulkk2){
       entry[0] = modulkk;
       entry[1] = modulkk2;
     }else{
       entry[0] = modulkk2;
       entry[1] = modulkk;
     }
     entry[2] = alpha;
     superSlabTF.push_back(entry);
   } 

   double areaPrimSub = calcArea(surfSub);
   double tempMCIA = nSub*areaPrimSub;
   vector<double> entry(14);
   // Compare bost list
   for(int i=0; i<superSlabSub.size(); i++){ 
     for(int j=0; j<superSlabTF.size(); j++){ 
       double errA = 100*(superSlabSub[i][0]-superSlabTF[j][0])/superSlabTF[j][0];
       double errB = 100*(superSlabSub[i][1]-superSlabTF[j][1])/superSlabTF[j][1];
       double errAlpha = 100*(superSlabSub[i][2]-superSlabTF[j][2])/superSlabTF[j][2]; 
      //cout << "SurperSub " << superSlabSub[i][0] << " " << superSlabSub[i][1] << " " << superSlabSub[i][2] << std::endl; 
      //cout << "SurperTF " << superSlabTF[i][0] << " " << superSlabTF[i][1] << " " << superSlabTF[i][2] << std::endl; 
	  //cout << " Errors " << errA << " " << errB << " " << errAlpha << std::endl;
      if(abs(errA) < 5.0 && abs(errB) < 5.0 && abs(errAlpha) < 5.0){ // IMPORTANT THRESHOLD !!!!!!!!
         entry[0]=tempMCIA;
         entry[1]=nSub;
         entry[2]=nTF;
         entry[3]=errA;
         entry[4]=errB;
         entry[5]=errAlpha;
         entry[6]=genMatSub[i][0][0];
         entry[7]=genMatSub[i][0][1];
         entry[8]=genMatSub[i][1][0];
         entry[9]=genMatSub[i][1][1];
         entry[10]=genMatTF[j][0][0];
         entry[11]=genMatTF[j][0][1];
         entry[12]=genMatTF[j][1][0];
         entry[13]=genMatTF[j][1][1];
         listR.push_back(entry);
       }
     }
   }
 }

 vector<vector<vector<int> > >matGen(int N){

   vector<vector<vector<int> > > genMatkk;

   for(int ii=1; ii<=N;ii++){
     if(N % ii == 0){
       int jj = N/ii;
       // i = ii and m == jj 
       for(int j=0; j<jj; j++){
         vector<vector<int> > kk(2, vector<int>(2));
         kk[0][0]=ii; 
         kk[0][1]=j; 
         kk[1][0]=0; 
         kk[1][1]=jj;
         genMatkk.push_back(kk); 
       }
     } 
   }   
   return genMatkk;
 }
 
 void readPlanes(){

   std::ifstream file_TF, file_Sub;
   file_TF.open("thinfilm.dat");
   file_Sub.open("substrate.dat");

   vector<int>  kkvec(3);
   while(file_TF >>  kkvec[0] >>  kkvec[1] >>  kkvec[2]){
     planesTF.push_back(kkvec);
   }
   while(file_Sub >>  kkvec[0] >>  kkvec[1] >>  kkvec[2]){
     planesSub.push_back(kkvec);
   }

 } 

 void printList(){
   //cout << "En print List " << std::endl;
   ofstream mylist;
   mylist.open ("list.dat");
   mylist << "[INTERFACES]***************************************************************************************" << std::endl;
   mylist << "[MATCHING INTERFACES]START" << std::endl;
   mylist << "    MCIA  N_TF  N_Sub   Err_a    Err_b  Err_alpha   MAT_TF     MAT_SUB      " << std::endl;
   for(int i=0; i<listR.size(); i++){
     mylist << setw(8) <<listR[i][0] << setw(6) <<  int(listR[i][1]) << setw(6) <<  int(listR[i][2])  <<  setprecision(4) << setw(10) <<  listR[i][3] << setprecision(4) <<  setw(10) << listR[i][4] << setprecision(3) << setw(10) <<  listR[i][5] << setw(6) <<  int(listR[i][6]) << setw(3) <<  int(listR[i][7]) << setw(3) << int(listR[i][8]) << setw(3) << int(listR[i][9]) << setw(6) <<  int(listR[i][10]) << setw(3) <<  int(listR[i][11]) << setw(3) <<  int(listR[i][12]) << setw(3) <<  int(listR[i][13]) << "  " << std::endl;
   } 
   mylist << "[MATCHING INTERFACES]STOP" << std::endl;
   mylist << "[INTERFACES]***************************************************************************************" << std::endl;
   mylist.close();
 }  

 void sortList(){
   //std::sort(listR.begin(), listR.end(),[](const std::vector<double>& a, const std::vector<double>& b){return a[0] < b[0];});
 }

 int main () {

   // Read thin film bulk POSCAR   
   
   // Read thin substrate bulk  POSCAR   

   // READ Planes
      //readPlanes;
   // Use Aflow for cutting planes

   // Read unit surface cell for thin film POSCAR 
   surfTF=readSurf("thinF.POSCAR");
   //cout << surfTF[0][0] << " " << surfTF[0][1] << " " << surfTF[1][0] << " "  << surfTF[1][1] << std::endl;
   // Calculate  unit surface cell  area for thin film POSCAR 
   double areaPrimTF = calcArea(surfTF);
   //cout << areaPrimTF << std::endl;
    // Read unit surface cell for substrate film POSCAR 
   surfSub=readSurf("subs.POSCAR");
   //cout << surfSub[0][0] << " " << surfSub[0][1] << " " << surfSub[1][0] << " "  << surfSub[1][1] << std::endl;
   // Calculate  unit surface cell  area for substrate POSCAR 
   double areaPrimSub = calcArea(surfSub);
   //cout << areaPrimSub << std::endl;

   // Check tentative MCIA and buld N1/N2 lookup table
   vector<vector<int> > listAreas;
   listAreas = calcAreasList(areaPrimTF,areaPrimSub);
   for(int i=0;i<listAreas.size();i++){
	cout << listAreas[i][0] << " " << listAreas[i][1] << std::endl;
     //Compare each potential matching areas
     compareSupperLattice(listAreas[i][0],listAreas[i][1]);
   }
   // Sort and print mathches
   sortList();
   if(listR.size()>0){
     printList();
   }
   
   return 0;
 }


