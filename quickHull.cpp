#ifndef QUICKHULL
#define QUICKHULL

#include <iostream>
#include <vector>
#include <string>
#include <utility> // for the std::pair
#include <sstream>
#include <fstream>
#include <cctype>
#include <math.h>
#include <cmath>
#include <set>
#include <unordered_set>

//https://stackoverflow.com/questions/8833938/is-the-stdset-iteration-order-always-ascending-according-to-the-c-specificat
// for knowledge on comparison structs

//struct for angle set comparison
struct cmpAngleStruct {
    bool operator() (std::pair<std::pair<double,double>,int> const & lhs, std::pair<std::pair<double,double>,int> const & rhs) const
    {
        if (lhs.first.first != rhs.first.first) {
        return lhs.first.first < rhs.first.first;
        // should be this: if (lhs.first.first == rhs.first.first)
        } else {
            return lhs.first.second < rhs.first.second;
        }
    }
};
// checks for y first, then x
struct cmpCVStruct {
    bool operator() (std::pair<std::pair<double,double>,int> const & lhs, std::pair<std::pair<double,double>,int> const & rhs) const
    {
        if (lhs.first.second != rhs.first.second) {
        return lhs.first.second < rhs.first.second;
        // should be this: if (lhs.first.second == rhs.first.second) (y1 == y2)
        } else {
            return lhs.first.first < rhs.first.first;
        }
    }
};


class QuickHull {
    private:
    std::vector<std::pair<double,double>> pointDex;
    std::set<std::pair<std::pair<double,double>,int>, cmpCVStruct> convexHull;
    std::set<std::pair<std::pair<double,double>,int>, cmpAngleStruct> angleOrder; // (angle, distance from center), input index
    int size;
    
    public:
    // Will add all points into the point vector, from the default txt file (points_100.txt)
    QuickHull(std::vector<std::pair<double,double>> inputs) {
        inputPoints(inputs);
        this->size = inputs.size();
    }
    // Adds all points in the specified txt file into the point vector.
    QuickHull(std::string pointsTxt){
        inputPoints(pointsTxt);
    }
    void inputPoints(std::string pointsTxt) {
        double xPoint;
        double yPoint;

        //file read-in credit given to cplusplus.com, at
        // https://www.cplusplus.com/doc/tutorial/files/
        std::string line;
        std::string segment;
        //is basically hobbies.txt
        std::ifstream myfile(pointsTxt);

        if (myfile.is_open()) {
            //initial run for size of thingy
            std::getline(myfile, line);
            this->size = std::stod(line);
            while ( std::getline(myfile, line) ) {
                // the line is now gotten, so time to push into the vector
                // inspiration from stackoverflow very nice, link is:
                // https://stackoverflow.com/questions/10058606/splitting-a-string-by-a-character
                std::stringstream doubleString(line);
                while (std::getline(doubleString, segment, ',')) {
                    // 2 states: a double with a space in front --> y, and a double without one --> x
                    if (std::isspace(segment[0])) {
                        yPoint = std::stod(segment.substr(1));
                    } else {
                        xPoint = std::stod(segment);
                    }
                }
                // add a point in
                addPoint(xPoint, yPoint);
            }
            myfile.close();
        } else {
            std::cerr << "Hobbies file could not be opened, check if the file exists or if program has necessary perms to access file.\n";
            exit(1);
        }

    }
    void inputPoints(std::vector<std::pair<double,double>> inputs) {
        for (int i = 0; i < (int)inputs.size(); i++) {
            addPoint(inputs[i].first, inputs[i].second);
        }
    }

    void addPoint(double xLoc, double yLoc) {
        std::pair<double, double> coordinate;
        coordinate.first = xLoc;
        coordinate.second = yLoc;
        pointDex.push_back(coordinate);
    }

    int whichSide(std::pair<double,double> lineP1, std::pair<double,double> lineP2, std::pair<double,double> point) {
        double val = (lineP1.second - point.second) * (lineP2.first - lineP1.first) -
              (lineP2.second - lineP1.second) * (lineP1.first - point.first);
  
        if (val > 0)
            return 1;
        if (val < 0)
            return -1;
        return 0;
    }
    // quickhull has the initial dividing step and the recursive steps, where we check the furthest point, make lines to that point,
    // and then do it again. 
    //On second thought, we aren't going to add another input for ignoreVals. It doesn't work with 1 million points.
    void quickHullMain(int size, std::pair<double,double> lineP1, std::pair<double, double> lineP2, int side, int indexP1, int indexP2, int counter) {
        int indx = -1; //for the side stuff
        double maxDist = 0; //for finding the furthest point
        double currDist;
        //I dont have the time to figure out wtf is going on with the stackoverflow
        if (counter > size) {
            return;
        } else {
            counter++;
        }
        // copy vals into new unordered set
        //std::unordered_set<int> nextIgnore (ignoreVals);


        // will check all points on whether its the right side, and/or whether its the furthest point away from the line so far.
        for (int i = 0; i < size; i++) {
            currDist = pointToLineDist(lineP1, lineP2, pointDex[i]);
            //hopefully this part reducess the overhead. My previous idea definitely did not lmao
            if (whichSide(lineP1, lineP2, pointDex[i]) == side && currDist > maxDist) {
                indx = i;
                maxDist = currDist;
            }
        }
        //last ditch force:

        convexHull.insert(std::make_pair(lineP2,indexP2));

        // ie, the thingy reaches the end w/o finding anything, so it puts the stuff into the convex hull set.
        // This is for redundancy.
        if (indx == -1) {
            convexHull.insert(std::make_pair(lineP1,indexP1));
            convexHull.insert(std::make_pair(lineP2,indexP2));
            return;
        }
        
        // check for debugging issues I'm having. Hopefully this will force things to close.
        if (pointDex[indx] == lineP1 || pointDex[indx] == lineP2) {
            convexHull.insert(std::make_pair(lineP1,indexP1));
            convexHull.insert(std::make_pair(lineP2,indexP2));
            return;
        }

        //recursion part
        // we need to check whether there are any points that need to be added on the left and right sides
        // of the point just added, hence this step:
        /*
        quickHullMain(size, lineP1, pointDex[indx], -whichSide(pointDex[indx], lineP1, lineP2), indexP1, indx);
        quickHullMain(size, lineP2, pointDex[indx], -whichSide(pointDex[indx], lineP2, lineP1), indexP2, indx);
        */
        quickHullMain(size, lineP1, pointDex[indx], side, indexP1, indx, counter);
        quickHullMain(size, lineP2, pointDex[indx], side, indexP2, indx, counter);
    }
    
    // really we only need how many points to use to construct the Convex Hull
    void quickHullHandler() {
        //first, wipe the stuff in convexHull:
        convexHull.clear();

        // 2 point Convex Hull not possible
        if (size < 3) {
            return; 
        }

        // find the first points, min and max x.
        std::pair<double,double> minX = pointDex[0];
        std::pair<double,double> maxX = pointDex[0];
        int minXIndx = 0;
        int maxXIndx = 0;

        for (int i = 0; i < size; i++) {
            if (minX > pointDex[i]) {
                minX = pointDex[i];
                minXIndx = i;
            } else if (minX == pointDex[i]) {
                if (minX.second > pointDex[i].second) {
                    minX = pointDex[i];
                    minXIndx = i;
                }
            } else if (maxX < pointDex[i]) {
                maxX = pointDex[i];
                maxXIndx = i;
            } else if (maxX == pointDex[i]) {
                if (maxX.second < pointDex[i].second) {
                    maxX = pointDex[i];
                    maxXIndx = i;
                }
            }
        }
        // one side
        quickHullMain(size, minX, maxX, 1, minXIndx, maxXIndx, 0);
    
        // and the other side
        quickHullMain(size, minX, maxX, -1, minXIndx, maxXIndx, 0);




    }

    //couts the points of the convex hull, counterclockwise order from the 0th vertex
    void outputCV() {
        //first convert all points into their polar coordinates:
        for (auto itr : convexHull) {
            pointToPolar(itr.first.first, itr.first.second, itr.second, convexHull.begin()->first);
        }
        // now everything is ordered in the set, based on angle and distance, so time to cout everything
        //cout number of stuff first
        std::cout << angleOrder.size() << std::endl;
        
        // and now output everything
        for (auto itr : angleOrder) {
        std::cout << itr.second << "," << pointDex[itr.second].first << "," << pointDex[itr.second].second << std::endl;
        }
    }

    void pointToPolar(double xLoc, double yLoc, int index, std::pair<double,double> center) {
        std::pair<double, double> polarCoords; //This should go (distance, angle) (angle from the 0th vertex)
        
        polarCoords.second = std::sqrt(std::pow((yLoc - center.second),2) + std::pow((xLoc - center.first),2)); 
        polarCoords.first = findAngle(xLoc, yLoc, center);
        angleOrder.insert(std::make_pair(polarCoords,index));
        //return polarCoords;
    }
    // This is the angle from the x axis to the point
    double findAngle(double xCoord, double yCoord, std::pair<double,double> center) {
        if ((xCoord - center.first) == 0) {
            if ((yCoord - center.second) > 0.0) {
                return 90.0;
            } else if (xCoord - center.first < 0.0){
                return 270.0;
            } else return 0.0; // really only for the 0th vertex for the CV for proper ordering
        }
        double unchangedAngle = std::atan((yCoord - center.second)/(xCoord - center.first));
        if ((xCoord - center.first) < 0) { //2nd quadrant (should also work for 3rd quarant)
            return (unchangedAngle + 180);
        } else if ((yCoord - center.second) < 0) { //4th quadrant
            return (unchangedAngle + 360);
        } else { // 1st quadrant
            return unchangedAngle;
        }
    }

    double pointToLineDist(std::pair<double,double> linePt1, std::pair<double,double> linePt2, std::pair<double,double> point) {
        double dist1 = std::abs((point.second - linePt1.second) * (linePt2.first - linePt1.first) -
                                (linePt2.second - linePt1.second) * (point.first - linePt1.first));
        double divide = std::sqrt(std::pow((linePt2.first -linePt1.first),2) + std::pow((linePt2.second -linePt1.second),2));
        double dist = dist1/divide;
        if (divide == 0) {
            return 0;
        }
        return dist;
    }

};

#endif