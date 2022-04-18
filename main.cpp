#include <iostream>
#include "quickHull.cpp"
int main() {
    //initial cin for the number of points:
    int size;
    std::cin >> size;
    std::string holdString = "";
    std::string segment;
    int isX = 0;
    double x;
    double y;
    std::vector<std::pair<double,double>> pointList;
    //cin will run from 0 to
    for (int i = 0; i < size; i++) {
        std::cin >> holdString;
        // now parse the string properly:
        std::stringstream doubleString(holdString);
                while (std::getline(doubleString, segment, ',')) {
                    // 2 states: a double with a space in front --> y, and a double without one --> x
                    if (isX == 0) {
                        isX++;
                        x = std::stod(segment);
                    } else if (isX == 1) {
                        isX++;
                        y = std::stod(segment);
                    } else if (isX > 1) {
                        std::cout << "bad input, you have more than 2 inputs, accepting first two inputs.\n";
                        break;
                    }
                }
                if (isX < 2) {
                    std::cout << "error, input has only 0 or 1 characters.\n";
                } else {
                isX = 0;
                // now set x and y.
                std::pair<double,double> point;
                point.first = x; point.second = y;
                //and add this to the vector
                pointList.push_back(point);
                }

    }
    // and now enter it into the quickHull object
    QuickHull hull(pointList);
    // then do the quickhull operation
    hull.quickHullHandler(size);
    // and now, the magic
    hull.outputCV();

    return 0;
}