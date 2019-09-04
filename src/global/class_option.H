/*!
  @file class_option.H
  @brief Declaration of a class/struct that holds options
  @author Robert Marskar
  @date   Aug. 2019
*/

#ifndef _CLASS_OPTION_
#define _CLASS_OPTION_

template <class T> class class_option {
public:

  class_option(){
  }
  
  class_option(T* a_value, const std::string a_prefix, const std::string a_description) {
    define(a_value, a_prefix, a_description);
  }
  
  ~class_option(){
  }

  void define(T* a_value, const std::string a_prefix, const std::string a_description){
    m_value       = a_value;
    m_prefix      = a_prefix;
    m_description = a_description;
  }

  // Get functions
  T* get_option()                {return m_value;};
  std::string get_prefix()      {return m_prefix;}
  std::string get_description() {return m_description;}

  // Set functions
  void set_value(T* a_value)                      {m_value = a_value;}
  void set_prefix(const std::string a_prefix)          {m_prefix = a_prefix;}
  void set_description(const std::string a_description){m_description = a_description;}

protected:
  
  T* m_value;
  std::string m_prefix;
  std::string m_description;
  
};
#endif
