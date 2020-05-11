#include "Logger.hpp"

namespace logging {

Logger * Logger::_p_logger = NULL;

Logger::Logger()
    : _verbosity(Debug)
{}

Logger::~Logger()
{
  if (_fout.is_open())
    _fout.close();
}

Logger & Logger::ref()
{
  if (!_p_logger)  // create a logger
  {
    static Logger logger_instance;
    _p_logger = &logger_instance;
  }
  return *_p_logger;
}

void Logger::set_file(const std::string file_name)
{
  if (_fout.is_open())
    _fout.close();
  _fout.open(file_name.c_str());
}

void Logger::set_verbosity(const LogLevel verbosity)
{
  _verbosity = verbosity;
}

void Logger::set_message_level(const LogLevel level)
{
  _msg_level = level;
}

void Logger::print_prefix_(const LogLevel level)
{
  // if (_fount.is_open())  // timestamp

  switch (level)
  {
    case Debug:
      _fout << "Debug   : ";
      break;
    case Message:
      _fout << "Message : ";
      break;
    case Brief:
      _fout << "Brief   : ";
      break;
    case Warning:
      std::cout << colors::BrightRed << "Warning: ";
      _fout << "Warning : ";
      break;
    case Critical:
      std::cout << colors::BrightRed << "Error: ";
      _fout << "Critical: ";
      break;
  }
}

}  // end namespace logger