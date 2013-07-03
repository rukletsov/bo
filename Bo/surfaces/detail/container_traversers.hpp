
/******************************************************************************

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
#include <memory>
#include <boost/shared_ptr.hpp>

namespace bo {
namespace surfaces {
namespace detail {

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

    TraverseRule(ContainerConstPtr contour, std::size_t index)
        : contour_(contour), current_index_(index)
    { }

    virtual void add(std::size_t offset) = 0;
    virtual Reference dereference() const = 0;
    virtual bool check_validity() const = 0;
    virtual SelfType* clone() const = 0;
    virtual std::size_t index()
    { return current_index_; }

    virtual ~TraverseRule()
    { }

    // Use this macro to define an implementation of the clone() function in descendants.
    #define BO_TRAVERSE_RULE_CLONE_IMPL \
        virtual SelfType* clone() const \
        { return new SelfType(*this); }

protected:
    ContainerConstPtr contour_;
    std::size_t current_index_;
};

// This traverse rule iterates the container forward from begin() to end().
template <typename Container>
class FwdOnePassTraverseRule: public TraverseRule<Container>
{
public:
    typedef FwdOnePassTraverseRule<Container> SelfType;
    typedef typename Container::const_iterator FwdIterator;
    typedef typename TraverseRule<Container>::Reference Reference;

    FwdOnePassTraverseRule(typename TraverseRule<Container>::ContainerConstPtr contour)
        : TraverseRule<Container>(contour, 0), end_it_(contour->end())
    { fwd_it_ = this->contour_->begin(); }

    virtual void add(std::size_t offset)
    { fwd_it_ += offset; this->current_index_ += offset; }

    virtual Reference dereference() const
    { return *fwd_it_; }

    virtual bool check_validity() const
    { return (fwd_it_ != end_it_); }

    virtual ~FwdOnePassTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    FwdIterator fwd_it_;
    FwdIterator end_it_;
};

// This traverse rule iterates the container backward from rbegin() to rend().
template <typename Container>
class BwdOnePassTraverseRule: public TraverseRule<Container>
{
public:
    typedef BwdOnePassTraverseRule<Container> SelfType;
    typedef typename Container::const_reverse_iterator BwdIterator;
    typedef typename TraverseRule<Container>::Reference Reference;

    BwdOnePassTraverseRule(typename TraverseRule<Container>::ContainerConstPtr contour)
        : TraverseRule<Container>(contour, contour->size() - 1), end_it_(contour->rend())
    { bwd_it_ = this->contour_->rbegin(); }

    virtual void add(std::size_t offset)
    { bwd_it_ += offset; this->current_index_ -= offset; }

    virtual Reference dereference() const
    { return *bwd_it_; }

    virtual bool check_validity() const
    { return (bwd_it_ != end_it_); }

    virtual ~BwdOnePassTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    BwdIterator bwd_it_;
    BwdIterator end_it_;
};

// This traverse rule iterates the container forward from the given item to back()
// and then from front() to the given item included.
template <typename Container>
class FwdCircuitTraverseRule: public FwdOnePassTraverseRule<Container>
{
public:
    typedef FwdCircuitTraverseRule<Container> SelfType;
    typedef typename TraverseRule<Container>::ContainerConstPtr ContainerConstPtr;

    FwdCircuitTraverseRule(ContainerConstPtr contour, std::size_t start_idx):
        FwdOnePassTraverseRule<Container>(contour), left_items_(contour->size())
    {
        this->fwd_it_ += start_idx;
        this->current_index_ = start_idx;
    }

    virtual void add(std::size_t offset)
    {
        // Do nothing for already invalidated traverser.
        if (this->check_validity())
        {
            // Check boundary condition.
            left_items_ -= offset;
            if (left_items_ >= 0)
            {
                // Check if we have to jump to the beginning of the container.
                std::size_t dist_to_end = this->end_it_ - this->fwd_it_;
                if (dist_to_end > offset)
                {
                    // Iterator's new position is before end().
                    this->fwd_it_ += offset;
                    this->current_index_ += offset;
                }
                else
                {
                    // Start iterating from the beginning, mind skipped items.
                    this->fwd_it_ = this->contour_->begin();
                    std::size_t additional_shift = offset - dist_to_end;
                    this->fwd_it_ += additional_shift;
                    this->current_index_ = additional_shift;
                }
            }
            else
            {
                // Invalidate traverser.
                this->fwd_it_ = this->end_it_;
            }
        }
    }

    virtual ~FwdCircuitTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    std::ptrdiff_t left_items_;
};

// This traverse rule iterates the container backwards from the given item to front()
// and then from back() to the given item included.
template <typename Container>
class BwdCircuitTraverseRule: public BwdOnePassTraverseRule<Container>
{
public:
    typedef BwdCircuitTraverseRule<Container> SelfType;
    typedef typename TraverseRule<Container>::ContainerConstPtr ContainerConstPtr;

    BwdCircuitTraverseRule(ContainerConstPtr contour, std::size_t start_idx):
        BwdOnePassTraverseRule<Container>(contour), left_items_(contour->size())
    {
        this->bwd_it_ += (this->contour_->size() - 1 - start_idx);
        this->current_index_ = start_idx;
    }

    virtual void add(std::size_t offset)
    {
        // Do nothing for already invalidated traverser.
        if (this->check_validity())
        {
            // Check boundary condition.
            left_items_ -= offset;
            if (left_items_ >= 0)
            {
                // Check if we have to jump to the end of the container.
                std::size_t dist_to_end = this->end_it_ - this->bwd_it_;
                if (dist_to_end > offset)
                {
                    // Iterator's new position is after rend().
                    this->bwd_it_ += offset;
                    this->current_index_ -= offset;
                }
                else
                {
                    // Start iterating from the end backwards, mind skipped itemd.
                    this->bwd_it_ = this->contour_->rbegin();
                    this->bwd_it_ += (offset - dist_to_end);
                    this->current_index_ = this->contour_->size() - 1 - offset + dist_to_end;
                }
            }
            else
            {
                // Invalidate traverser.
                this->bwd_it_ = this->end_it_;
            }
        }
    }

    virtual ~BwdCircuitTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    std::ptrdiff_t left_items_;
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
    typedef std::auto_ptr<TraverseRuleType> TraverseRulePtr;

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

    std::size_t index() const
    {
        return rule_->index();
    }

    // Necessary for assignment operator.
    void swap(SelfType& other)
    {
        std::swap(valid_, other.valid_);
        std::swap(rule_, other.rule_);
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
    typedef std::auto_ptr<TraverseRuleType> TraverseRulePtr;

    #define BO_RULE_UNARY_FACTORY_FUNCTION(RuleName)                        \
        static TraverseRulePtr RuleName(ContainerConstPtr container_ptr)    \
        {                                                                   \
            typedef RuleName##TraverseRule<Container> Rule;                 \
            TraverseRulePtr rule_ptr(new Rule(container_ptr));              \
            return rule_ptr;                                                \
        }

    #define BO_RULE_BINARY_FACTORY_FUNCTION(RuleName)                       \
        static TraverseRulePtr RuleName(ContainerConstPtr container_ptr,    \
                                        std::size_t start_idx)              \
        {                                                                   \
            typedef RuleName##TraverseRule<Container> Rule;                 \
            TraverseRulePtr rule_ptr(new Rule(container_ptr, start_idx));   \
            return rule_ptr;                                                \
        }

    BO_RULE_UNARY_FACTORY_FUNCTION(FwdOnePass)
    BO_RULE_UNARY_FACTORY_FUNCTION(BwdOnePass)
    BO_RULE_BINARY_FACTORY_FUNCTION(FwdCircuit)
    BO_RULE_BINARY_FACTORY_FUNCTION(BwdCircuit)

    // TODO: create enums and true Create dispatcher.

    static TraverseRulePtr Create(ContainerConstPtr container_ptr, bool is_forward)
    {
        TraverseRulePtr retvalue = (is_forward ? FwdOnePass(container_ptr) :
                                                 BwdOnePass(container_ptr));
        return retvalue;
    }

    static TraverseRulePtr Create(ContainerConstPtr container_ptr,
                                std::size_t start_idx, bool is_forward)
    {
        TraverseRulePtr retvalue = (is_forward ? FwdCircuit(container_ptr, start_idx) :
                                                 BwdCircuit(container_ptr, start_idx));
        return retvalue;
    }
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // CONTAINER_TRAVERSERS_HPP_79A9B7BE_23EA_47D8_AAB9_2019FE8B4EF3
