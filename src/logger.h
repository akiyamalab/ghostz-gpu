/*
 * logger.h
 *
 *  Created on: 2011/04/11
 *      Author: shu
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <iostream>
#include <ostream>
#include <string>

class Logger {
public:

  static Logger* GetInstance() {
    static Logger instance;
    return &instance;
 }

  void ErrorLog(std::string message) {
//#pragma omp critical(lock_)
    std::cout << "error : " << message << std::endl;
  }

  void WarningLog(std::string message) {
//#pragma omp critical(lock_)
    std::cout << "warning : " << message << std::endl;
  }

  void Log(std::string message) {
//#pragma omp critical(lock_)
    std::cout << message << std::endl;
  }

private:
  Logger()
  {
  }
  ~Logger() {
  }
  Logger(const Logger& rhs);
  Logger operator=(const Logger& rhs);
};

#endif /* LOGGER_H_ */
