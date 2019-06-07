//
// Created by Yue Hao on 7/15/16.
//

#ifndef EPIC_QUAD_TREE_H
#define EPIC_QUAD_TREE_H
#include <vector>
#include <cmath>
#include <complex>
const unsigned long MAX_TREE_CAPACITY=10;

class CPoint{
public:
    std::complex<double> coordinate;

    CPoint():coordinate(0.0,0.0) {}
    CPoint(const double& xx, const double& yy): coordinate(xx, yy){}
    CPoint(const std::complex<double> coord) : coordinate(coord){}
    const double x() const {return coordinate.real();}
    const double y() const {return coordinate.imag();}
};

class box{
public:
    CPoint center;
    double hw;  // Half width of this box

    box(const CPoint& c, const double& width):center(c),hw(width){}
    box sub_box(const int& index){
        double cx,cy;
        if (index==1){cx=-0.5; cy=0.5;}
        else if (index==2){cx=+0.5; cy=0.5;}
        else if (index==3){cx=-0.5; cy=-0.5;}
        else {cx=0.5; cy=-0.5;}
        double new_x=center.x()+cx*hw;
        double new_y=center.y()+cy*hw;
        return box(CPoint(new_x,new_y),hw/2.0);
    }

    bool in_box(const CPoint& p){
        if (p.x()<=center.x()+hw && p.x()>center.x()-hw && p.y()<=center.y()+hw && p.y()>center.y()-hw) return true;
        return false;
    }
    bool in_box(const double& x, const double & y){
        return in_box(CPoint(x,y));
    }
    bool neighbor(box other){
        double disx=fabs(center.x()-other.center.x());
        double disy=fabs(center.y()-other.center.y());
        if (disx>hw+other.hw || disy>hw+other.hw) return false;
        return true;
    }


};
class CQuad_tree {
private:
    const unsigned long MAX_CAPACITY=10;

    box boundary;
    unsigned level;


    std::vector<CPoint> source;
    std::vector<CPoint> target;
    std::vector<CQuad_tree *> interaction_list;
    CQuad_tree * parent;
    CQuad_tree * q1; //nw
    CQuad_tree * q2; //ne
    CQuad_tree * q3; //sw
    CQuad_tree * q4; //se
public:
    static unsigned max_index;
    static std::vector<CQuad_tree *>  pointer_list;
    unsigned index;
    CQuad_tree(const box& boundary, const unsigned& level=0):
            boundary(boundary), MAX_CAPACITY(MAX_TREE_CAPACITY), level(level){
        index=max_index;
        max_index+=1;
        pointer_list.push_back(this);
        this->parent = nullptr;
        this->q1 = nullptr;
        this->q2 = nullptr;
        this->q3 = nullptr;
        this->q4 = nullptr;
        this->interaction_list.reserve(27);
        this->source.reserve(MAX_CAPACITY);
    }
    ~CQuad_tree(){
        delete q1;
        delete q2;
        delete q3;
        delete q4;
    }
    void subdivide();
    bool insert(CPoint p);
    void create_interaction_list();
    std::vector<double> get_boundary();

};
unsigned CQuad_tree::max_index = 0;
std::vector<CQuad_tree *> CQuad_tree::pointer_list={};

#endif //EPIC_QUAD_TREE_H
