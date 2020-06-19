#pragma once 

#include <vector>
#include <CGAL/Random.h>

struct directions{
    double _x, _y, _z;
};


class RaysGenerate{
private:
public:
    RaysGenerate(int n);
    RaysGenerate(int n, unsigned int seed);
    ~RaysGenerate();

    int _n;
    unsigned int _seed;
    std::vector<directions> rayDirections;
    std::vector<directions> normalisedRayDirections;
};

RaysGenerate::RaysGenerate(int n): _n(n){

    CGAL::Random rand;

    for(size_t i=0; i!=n; ++i){

        directions direction;
        direction._x = rand.get_double(-1.0, 1.0);
        direction._y = rand.get_double(-1.0, 1.0); 
        direction._z = rand.get_double(-1.0, 1.0); 
        rayDirections.push_back(direction);

        directions normalisedDirection;
        normalisedDirection._x = direction._x/ sqrt(pow(direction._x, 2) + pow(direction._y, 2) + pow(direction._z, 2));
        normalisedDirection._y = direction._y/ sqrt(pow(direction._x, 2) + pow(direction._y, 2) + pow(direction._z, 2));
        normalisedDirection._z = direction._z/ sqrt(pow(direction._x, 2) + pow(direction._y, 2) + pow(direction._z, 2));
        normalisedRayDirections.push_back(normalisedDirection);
    }
}

RaysGenerate::RaysGenerate(int n, unsigned int seed): _n(n), _seed(seed){

    CGAL::Random rand(seed);

    for(size_t i=0; i!=n; ++i){

        directions direction;
        direction._x = rand.get_double(-1.0, 1.0);
        direction._y = rand.get_double(-1.0, 1.0); 
        direction._z = rand.get_double(-1.0, 1.0); 
        rayDirections.push_back(direction);

        directions normalisedDirection;
        normalisedDirection._x = direction._x/ sqrt(pow(direction._x, 2) + pow(direction._y, 2) + pow(direction._z, 2));
        normalisedDirection._y = direction._y/ sqrt(pow(direction._x, 2) + pow(direction._y, 2) + pow(direction._z, 2));
        normalisedDirection._z = direction._z/ sqrt(pow(direction._x, 2) + pow(direction._y, 2) + pow(direction._z, 2));
        normalisedRayDirections.push_back(normalisedDirection);

    }
}

RaysGenerate::~RaysGenerate(){}
