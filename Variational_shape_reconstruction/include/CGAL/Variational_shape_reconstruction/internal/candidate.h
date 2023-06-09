#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "types.h"
namespace qem {

template <class Handle>
class Candidate
{
private:
    Handle      m_handle;
    int         m_index;
    double      m_loss;

public:

	Candidate(  Handle handle,
                int    index,
                double loss):
                m_handle(handle),
                m_index(index),
                m_loss(loss)
    {}

    ~Candidate() {} 

public:

	const double& loss() const { return m_loss; }
    double& loss() { return m_loss; }

	const Handle handle() const { return m_handle; }
    Handle handle() { return m_handle; }

    const int& index() const { return m_index; }
    int& index() { return m_index; }

    bool operator==(const Candidate<Handle>& obj) const
    {
        if(m_handle == obj.handle() && m_index == obj.index())
            return true;
        else
            return false;
    }
};  // end of class Candidate

template <class Handle>
class PCandidate
{
private:
    Handle      m_handle;
    Point       m_optimal;
    double      m_loss;

public:

	PCandidate( Handle handle,
                Point  optimal,
                double loss):
                m_handle(handle),
                m_optimal(optimal),
                m_loss(loss)
    {}

    ~PCandidate() {} 

public:

	const double& loss() const { return m_loss; }
    double& loss() { return m_loss; }

	const Handle handle() const { return m_handle; }
    Handle handle() { return m_handle; }
	Point optimal() { return m_optimal; }

    bool operator==(const PCandidate<Handle>& obj) const
    {
        if(m_handle == obj.handle())
            return true;
        else
            return false;
    }
};  // end of class PCandidate


template <class CCandidate>
struct Candidate_more
{
    // c1 > c2 means c1 is less prioritised over c2
    bool operator()(const CCandidate& c1, const CCandidate& c2)
    {
        return c1.loss() > c2.loss();
    }

}; // end of struct Candidate_more

template <class Handle>
struct HashHandle
{
public:
    size_t operator()(const Handle& elem) const
    {
        return elem;
    }
}; // end of HashHandle

struct HashPairIndex
{
    std::size_t operator () (const std::pair<int, int> &p) const 
    {
        std::size_t hash = 0;
        boost::hash_combine(hash, std::min(p.first, p.second));
        boost::hash_combine(hash, std::max(p.first, p.second));
        return hash;
    }
}; // end of HashPairIndex

struct EqualPairIndex 
{
    bool operator()(const std::pair<int, int> &p1, const std::pair<int, int> &p2) const 
    {
        if((p1.first == p2.first) && (p1.second == p2.second))
            return true;

        if((p1.first == p2.second) && (p1.second == p2.first))
            return true;

        return false;
    }
};

struct HashIntSet
{
public:
    size_t operator()(const std::set<int>& candidate) const
    {
        std::size_t hash = 0;
        for(auto &elem: candidate)
            boost::hash_combine(hash, elem);

        return hash;
    }
}; // end of HashIntSet

} // namespace qem

#endif