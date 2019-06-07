//
// Created by Yue Hao on 7/15/16.
//

#include "quad_tree.h"

void CQuad_tree::subdivide() {
    q1=new CQuad_tree(this->boundary.sub_box(1), level=this->level+1);
    q1->parent=this;
    q2=new CQuad_tree(this->boundary.sub_box(2), level=this->level+1);
    q2->parent=this;
    q3=new CQuad_tree(this->boundary.sub_box(3), level=this->level+1);
    q3->parent=this;
    q4=new CQuad_tree(this->boundary.sub_box(4), level=this->level+1);
    q4->parent=this;
}
bool CQuad_tree::insert(CPoint p) {
    if (this->boundary.in_box(p)==false) return false;

    if (this->source.size() >= MAX_CAPACITY){
        this->subdivide();
        if (q1->insert(p)) return true;
        if (q2->insert(p)) return true;
        if (q3->insert(p)) return true;
        if (q4->insert(p)) return true;

    }
    else{
        this->source.push_back(p);
        return true;
    }
}

void CQuad_tree::create_interaction_list() {
    
}
std::vector<double> CQuad_tree::get_boundary() {
    double leftx = this->boundary.center.x()-this->boundary.hw;
    double rightx = this->boundary.center.x()+this->boundary.hw;
    double topy = this->boundary.center.x()+this->boundary.hw;
    double bottomy = this->boundary.center.x()-this->boundary.hw;
    return std::vector<double> {leftx, rightx, bottomy, topy};
}
