/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
          Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
          Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
          IMDEA Materiales Institute
 
**********************************************************************/
#ifndef segmenteddata_h
#define segmenteddata_h

typedef enum  {
    MINIMA,
    ONESADDLE,
    TWOSADDLE,
    MAXIMA
} criticalType;


/*
 
criticalPoint class stores the value of :
 - coordinates of the critical point in coordinates[3]
 - the type of critical point
 - and the scalar value at the critical point
 
 */


class criticalPoint {
    
public:
    criticalPoint() :
    ctype(MINIMA), scalarValue(0.0)
    {
        for (size_t i = 0;i<3;i++){coordinates[i] = 0.0;};
        
    };
    
    ~criticalPoint() {};
    
    void                                   setCoordinates(const double p[3]) {
                                                coordinates[0] = p[0];
                                                coordinates[1] = p[1];
                                                coordinates[2] = p[2];
    };
    
    void                                   getCoordinates(double p[3]) {
                                                p[0] = coordinates[0];
                                                p[1] = coordinates[1];
                                                p[2] = coordinates[2];
    };
        
    void                                    setCriticalType (size_t p){ ctype = static_cast<criticalType>(p);};
    criticalType                            getCriticalType (){return ctype;};
    
    void                                    setScalarValue (const double c) { scalarValue = c;};
    double                                  getScalarValue () {return scalarValue;};
    
    
private:
    
    double                                  coordinates[3];
    criticalType                            ctype;
    double                                  scalarValue;
    
};



/* To traverse the data, the innermost loop is of dimension nx, next on ny and then nz.
  i.e. (0,0,0), (1,0,0) ... (nx-1,0,0), (0,1,0) ... (nx-1, ny-1, nz-1) */

class segmentedData {
    
public:
    
    segmentedData() :                               coordinates(), ascendingManifoldID(), descendingManifoldID(),scalarValue() {};
    ~segmentedData() {};
    
    
    
    void                                            setDimensions(size_t nx, size_t ny, size_t nz){
                                                        coordinates.resize(nx*ny*nz, vector<double>(3));
                                                        ascendingManifoldID.resize(nx*ny*nz);
                                                        descendingManifoldID.resize(nx*ny*nz);
                                                        scalarValue.resize(nx*ny*nz);
        
                                                    for (size_t i = 0; i < nx*ny*nz; i++){
                                                        ascendingManifoldID[i] = 0;
                                                        descendingManifoldID[i] = 0;
                                                        scalarValue[i] = 0.0;
                                                        for (size_t j = 0; j < 3; j++){
                                                            coordinates[i][j] = 0.0;
                                                        }
                                                    }
        
    };
    
    void                                            setCoordinates(size_t i, const double p[3]) {
                                                        coordinates[i][0] = p[0];
                                                        coordinates[i][1] = p[1];
                                                        coordinates[i][2] = p[2];
    };
    
    void                                            getCoordinates(size_t i, double p[3]){
                                                        p[0] = coordinates[i][0];
                                                        p[1] = coordinates[i][1];
                                                        p[2] = coordinates[i][2];
    };
    
    void                                            setAscendingManifoldID(size_t i, size_t ID) {ascendingManifoldID[i] = ID;};
    size_t                                          getAscendingManifoldID(size_t i) {return ascendingManifoldID[i];};
    
    void                                            setDescendingManifoldID(size_t i, size_t ID) {descendingManifoldID[i] = ID;};
    size_t                                          getDescendingManifoldID(size_t i) {return descendingManifoldID[i];};
    
    void                                            setScalarValue(size_t i, double value) {scalarValue[i] = value;};
    double                                          getScalarValue(size_t i) {return scalarValue[i];};
    
    size_t                                          getNumberOfPoints() {return scalarValue.size();};
    
    
private:
    
    std::vector<std::vector<double>>                coordinates;
    std::vector<size_t>                             ascendingManifoldID;
    std::vector<size_t>                             descendingManifoldID;
    std::vector<double>                             scalarValue;
    
};


#endif /* segmenteddata_h */
