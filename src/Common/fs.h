#pragma once
#include <string>

bool file_exists(const std::string fname);
bool dir_exists(const std::string fname);
bool exists(const std::string fname);
bool file_can_read(const std::string fname);
bool file_can_create(const std::string fname);
bool dir_can_write(const std::string pname);
std::string to_absolute(const std::string xname);
std::string to_canonical(const std::string xname);
std::string dirname(const std::string xname);
bool create_directory(std::string base, std::string extend);
bool regular_file(const std::string fname);
bool starts_with(std::string fname, std::string pattern);
void rename(std::string lhs, std::string rhs);
void remove(std::string fname, bool recursive = false);
std::string replace_ext(std::string fname, std::string new_ext);
