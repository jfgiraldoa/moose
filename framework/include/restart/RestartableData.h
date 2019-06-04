//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "DataIO.h"

// C++ includes
#include <vector>

// Forward declarations
class RestartableDataValue;

/**
 * Abstract definition of a RestartableData value.
 */
class RestartableDataValue
{
public:
  /**
   * Constructor
   * @param name The full (unique) name for this piece of data.
   * @param context 'typeless' pointer to user-specific data.
   */
  RestartableDataValue(std::string name, void * context) : _name(name), _context(context) {}

  /**
   * Destructor.
   */
  virtual ~RestartableDataValue() = default;

  /**
   * String identifying the type of parameter stored.
   * Must be reimplemented in derived classes.
   */
  virtual std::string type() = 0;

  /**
   * The full (unique) name of this particular piece of data.
   */
  std::string name() { return _name; }

  /**
   * A context pointer for helping with load / store.
   */
  void * context() { return _context; }

  virtual void swap(RestartableDataValue * rhs) = 0;

  // save/restore in a file
  virtual void store(std::ostream & stream) = 0;
  virtual void load(std::istream & stream) = 0;

protected:
  /// The full (unique) name of this particular piece of data.
  std::string _name;

  /// A context pointer for helping with load and store
  void * _context;
};

/**
 * Concrete definition of a parameter value
 * for a specified type.
 */
template <typename T>
class RestartableData : public RestartableDataValue
{
public:
  /**
   * Constructor
   * @param name The full (unique) name for this piece of data.
   * @param context 'typeless' pointer to user-specific data.
   */
  RestartableData(std::string name, void * context) : RestartableDataValue(name, context)
  {
    _value_ptr = libmesh_make_unique<T>();
  }

  /**
   * @returns a read-only reference to the parameter value.
   */
  T & get() { return *_value_ptr; }

  /**
   * @returns a writable reference to the parameter value.
   */
  T & set() { return *_value_ptr; }

  /**
   * String identifying the type of parameter stored.
   */
  virtual std::string type() override;

  /**
   * Swap
   */
  virtual void swap(RestartableDataValue * rhs) override;

  /**
   * Store the RestartableData into a binary stream
   */
  virtual void store(std::ostream & stream) override;

  /**
   * Load the RestartableData from a binary stream
   */
  virtual void load(std::istream & stream) override;

private:
  /// Stored value.
  std::unique_ptr<T> _value_ptr;
};

// ------------------------------------------------------------
// RestartableData<> class inline methods
template <typename T>
inline std::string
RestartableData<T>::type()
{
  return typeid(T).name();
}

template <typename T>
inline void
RestartableData<T>::swap(RestartableDataValue * libmesh_dbg_var(rhs))
{
  mooseAssert(rhs, "Assigning NULL?");
  //  _value.swap(cast_ptr<RestartableData<T>*>(rhs)->_value);
}

template <typename T>
inline void
RestartableData<T>::store(std::ostream & stream)
{
  T & tmp = *_value_ptr;
  storeHelper(stream, tmp, _context);
}

template <typename T>
inline void
RestartableData<T>::load(std::istream & stream)
{
  loadHelper(stream, *_value_ptr, _context);
}

/**
 * Container for storing material properties
 */
using RestartableDatas = std::vector<std::map<std::string, std::unique_ptr<RestartableDataValue>>>;

