#ifndef __LOAD_HPP__
#define __LOAD_HPP__

#include <memory>
#include <armadillo>
#include "element.hpp"

/* TODO: write EVERYTHING more generally (Element instead of Q4, etc.) */
class Traction
{
public:
    Load(const Q4& elem, const arma::vec& traction, const const unsigned int e) 
        : elem_ptr(std::shared_ptr<Q4>(&elem)), 
        trac_ptr(std::shared_ptr<arma::vec>(&traction)), edge(e) {}
    
    std::shared_ptr<Q4> elem_ptr;
    std::shared_ptr<arma::vec> trac_ptr;
    unsigned int e;
};

#endif /* __LOAD_HPP__ */
