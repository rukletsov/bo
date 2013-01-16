
/******************************************************************************

  container_traversers.hpp, v 1.1.2 2013.01.16

  Various traversers for containers with contour points. Intended for use in
  triangulation algorithms.

  Copyright (c) 2013
  Alexander Rukletsov <rukletsov@gmail.com>
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1.  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
  2.  Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
  SUCH DAMAGE.

*******************************************************************************/

#ifndef CONTAINER_TRAVERSERS_HPP_79A9B7BE_23EA_47D8_AAB9_2019FE8B4EF3
#define CONTAINER_TRAVERSERS_HPP_79A9B7BE_23EA_47D8_AAB9_2019FE8B4EF3

#include <iterator>
#include <stdexcept>
#include <algorithm>
#include <boost/shared_ptr.hpp>

namespace bo {

// Base class for traverse rules. A traverse rule implements the particular traversing
// method for the points container, e.g. from the end backwards.
template <typename Container>
class TraverseRule
{
public:
    typedef TraverseRule<Container> SelfType;
    typedef boost::shared_ptr<const Container> ContainerConstPtr;
    typedef typename Container::const_iterator ContainerConstIterator;
    typedef typename std::iterator_traits<ContainerConstIterator>::reference Reference;

    TraverseRule(ContainerConstPtr contour): contour_(contour)
    { }

    virtual void add(std::size_t offset) = 0;
    virtual Reference dereference() const = 0;
    virtual bool check_validity() const = 0;
    virtual SelfType* clone() const = 0;

    virtual ~TraverseRule()
    { }

    // TODO: add more macros to extract common code from virtual functions.

    // Use this macro to define an implementation of the clone() function in descendants.
    #define BO_TRAVERSE_RULE_CLONE_IMPL \
        virtual SelfType* clone() const \
        { return new SelfType(*this); }

protected:
    ContainerConstPtr contour_;
};

// This traverse rule iterates the container forward from begin() to end().
template <typename Container>
class FwdOnePassTraverseRule: public TraverseRule<Container>
{
public:
    typedef FwdOnePassTraverseRule<Container> SelfType;
    typedef typename Container::const_iterator FwdIterator;

    FwdOnePassTraverseRule(ContainerConstPtr contour): TraverseRule(contour)
    { fwd_it_ = contour_->begin(); }

    virtual void add(std::size_t offset)
    { fwd_it_ += offset; }

    virtual Reference dereference() const
    { return *fwd_it_; }

    virtual bool check_validity() const
    { return (fwd_it_ != contour_->end()); }

    virtual ~FwdOnePassTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    FwdIterator fwd_it_;
};

// This traverse rule iterates the container backward from rbegin() to rend().
template <typename Container>
class BwdOnePassTraverseRule: public TraverseRule<Container>
{
public:
    typedef BwdOnePassTraverseRule<Container> SelfType;
    typedef typename Container::const_reverse_iterator BwdIterator;

    BwdOnePassTraverseRule(ContainerConstPtr contour): TraverseRule(contour)
    { bwd_it_ = contour_->rbegin(); }

    virtual void add(std::size_t offset)
    { bwd_it_ += offset; }

    virtual Reference dereference() const
    { return *bwd_it_; }

    virtual bool check_validity() const
    { return (bwd_it_ != contour_->rend()); }

    virtual ~BwdOnePassTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    BwdIterator bwd_it_;
};

// This traverse rule iterates the container forward from the given item to back()
// and then from front() to the given item included.
template <typename Container>
class FwdCircuitTraverseRule: public TraverseRule<Container>
{
public:
    typedef FwdCircuitTraverseRule<Container> SelfType;
    typedef typename Container::const_iterator FwdIterator;

    FwdCircuitTraverseRule(ContainerConstPtr contour, std::size_t start_idx):
        TraverseRule(contour), left_items_(contour->size())
    {
        fwd_it_ = contour_->begin();
        fwd_it_ += start_idx;
        valid_ = (fwd_it_ != contour_->end());
    }

    virtual void add(std::size_t offset)
    {
        left_items_ -= offset;
        if (left_items_ < 0)
        {
            valid_ = false;
            fwd_it_ = contour_->end();
            return;
        }

        std::ptrdiff_t dist_to_end = contour_->end() - fwd_it_;
        if (dist_to_end > 1)
            // can be progressed up to the last element, before end()
            fwd_it_ += offset;
        else
        {
            // move to the end and then start from begin()
            offset -= dist_to_end;
            fwd_it_ = contour_->begin();
            // progress iterator by remaining elements
            fwd_it_ += offset;
        }
    }

    virtual Reference dereference() const
    { return *fwd_it_; }

    virtual bool check_validity() const
    { return valid_; }

    virtual ~FwdCircuitTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    FwdIterator fwd_it_;
    std::ptrdiff_t left_items_;
    bool valid_;
};

// This traverse rule iterates the container backwards from the given item to front()
// and then from back() to the given item included.
template <typename Container>
class BwdCircuitTraverseRule: public TraverseRule<Container>
{
public:
    typedef BwdCircuitTraverseRule<Container> SelfType;
    typedef typename Container::const_reverse_iterator BwdIterator;

    BwdCircuitTraverseRule(ContainerConstPtr contour, std::size_t start_idx):
        TraverseRule(contour), left_items_(contour->size())
    {
        bwd_it_ = contour_->rbegin();
        bwd_it_ += (contour_->size() - 1 - start_idx);
        valid_ = (bwd_it_ != contour_->rend());
    }

    virtual void add(std::size_t offset)
    {
        left_items_ -= offset;
        if (left_items_ < 0)
        {
            valid_ = false;
            bwd_it_ = contour_->rend();
            return;
        }

        std::ptrdiff_t dist_to_end = contour_->rend() - bwd_it_;
        if (dist_to_end > 1)
            // can be progressed up to the first element, after rend()
            bwd_it_ += offset;
        else
        {
            // move to the front and then start from rbegin()
            offset -= dist_to_end;
            bwd_it_ = contour_->rbegin();
            // progress iterator by remaining elements
            bwd_it_ += offset;
        }
    }

    virtual Reference dereference() const
    { return *bwd_it_; }

    virtual bool check_validity() const
    { return valid_; }

    virtual ~BwdCircuitTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    BwdIterator bwd_it_;
    std::ptrdiff_t left_items_;
    bool valid_;
};

// This class represents a constant container traverser. The actual behaviour depends
// on the traverse rule specified. Turns into invalid state if the traverse rule
// indicates that no more elements can be accessed in the container.
template <typename Container>
class ContainerConstTraverser
{
public:
    typedef ContainerConstTraverser<Container> SelfType;

    typedef typename Container::const_iterator ContainerConstIterator;
    typedef typename std::iterator_traits<ContainerConstIterator>::reference Reference;

    typedef TraverseRule<Container> TraverseRuleType;
    typedef boost::shared_ptr<TraverseRuleType> TraverseRulePtr;

    ContainerConstTraverser(): valid_(false)
    { }

    ContainerConstTraverser(TraverseRulePtr rule_ptr): rule_(rule_ptr)
    {
        valid_ = rule_->check_validity();
    }

    // Custom copy c-tor for correct copying.
    ContainerConstTraverser(const SelfType& other): valid_(other.valid_),
        rule_(other.rule_->clone())
    { }

    // Custom assignment operator for correct copying.
    SelfType& operator=(SelfType other)
    {
        other.swap(*this);
        return *this;
    }

    SelfType& operator++()
    {
        return ((*this) += 1);
    }

    SelfType& operator+=(std::size_t offset)
    {
        if (is_valid())
        {
            rule_->add(offset);
            valid_ = rule_->check_validity();
        }

        return (*this);
    }

    SelfType operator+(std::size_t offset) const
    {
        SelfType temp = *this;
        return
            (temp += offset);
    }

    Reference operator*() const
    {
        if (!is_valid())
            throw std::logic_error("Cannot dereference invalid CommonConstIterator.");

        return
            rule_->dereference();
    }

    // Necessary for assignment operator.
    void swap(SelfType& other)
    {
        std::swap(valid_, other.valid_);
        rule_.swap(other.rule_);
    }

    bool is_valid() const
    {
        return valid_;
    }

private:
    bool valid_;
    TraverseRulePtr rule_;
};

// Factory for traverse rules. Simplifies construction of traversers with particular
// traverse rules.
template <typename Container>
struct TraverseRuleFactory
{
    typedef boost::shared_ptr<const Container> ContainerConstPtr;
    typedef TraverseRule<Container> TraverseRuleType;
    typedef boost::shared_ptr<TraverseRuleType> TraverseRulePtr;

    #define BO_RULE_FACTORY_FUNCTION(RuleName)                              \
        static TraverseRulePtr RuleName(ContainerConstPtr container_ptr)    \
        {                                                                   \
            typedef RuleName##TraverseRule<Container> Rule;                 \
            TraverseRulePtr rule_ptr(new Rule(container_ptr));              \
            return rule_ptr;                                                \
        }


    static TraverseRulePtr FwdCircuit(ContainerConstPtr container_ptr, std::size_t start_idx)
    {
        typedef FwdCircuitTraverseRule<Container> Rule;
        TraverseRulePtr rule_ptr(new Rule(container_ptr, start_idx));
        return rule_ptr;
    }

    static TraverseRulePtr BwdCircuit(ContainerConstPtr container_ptr, std::size_t start_idx)
    {
        typedef BwdCircuitTraverseRule<Container> Rule;
        TraverseRulePtr rule_ptr(new Rule(container_ptr, start_idx));
        return rule_ptr;
    }

    BO_RULE_FACTORY_FUNCTION(FwdOnePass)
    BO_RULE_FACTORY_FUNCTION(BwdOnePass)

    static TraverseRulePtr Create(ContainerConstPtr container_ptr, bool is_forward)
    {
        return (is_forward ? FwdOnePass(container_ptr) : BwdOnePass(container_ptr));
    }

    static TraverseRulePtr Create(ContainerConstPtr container_ptr,
                                std::size_t start_idx, bool is_forward)
    {
        if ((start_idx == 0) || (start_idx == container_ptr->size() - 1))
            return Create(container_ptr, is_forward);
        else
            return (is_forward ? FwdCircuit(container_ptr, start_idx) :
                                 BwdCircuit(container_ptr, start_idx));
    }
};

} // namespace bo

#endif // CONTAINER_TRAVERSERS_HPP_79A9B7BE_23EA_47D8_AAB9_2019FE8B4EF3
