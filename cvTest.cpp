#include <iostream>
#include "quickHull.cpp"

int main () {
    std::cout << "testing 20 point CV" << std::endl;

    std::string pointTxt;
    pointTxt = "prog1/points_20.txt";

    QuickHull hullz (pointTxt);
    hullz.quickHullHandler();
    hullz.outputCV();

    return 0;
}