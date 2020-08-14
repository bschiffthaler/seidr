/* * Seidr - Create and operate on gene crowd networks
 * Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
 *
 * This file is part of Seidr.
 *
 * Seidr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Seidr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include <string>

bool
file_exists(const std::string& fname);
bool
dir_exists(const std::string& fname);
bool
exists(const std::string& fname);
bool
file_can_read(const std::string& fname);
bool
file_can_create(const std::string& fname);
bool
dir_can_write(const std::string& pname);
std::string
to_absolute(const std::string& xname);
std::string
to_canonical(const std::string& xname);
std::string
dirname(const std::string& xname);
std::string
basename(const std::string& xname);
bool
create_directory(const std::string& path);
bool
create_directory(const std::string& base, const std::string& extend);
bool
regular_file(const std::string& fname);
bool
starts_with(const std::string& fname, const std::string& pattern);
void
rename(const std::string& lhs, const std::string& rhs);
void
remove(const std::string& fname, bool recursive = false);
std::string
replace_ext(const std::string& fname, const std::string& new_ext);
std::string
tempfile(const std::string& tempdir = "");
std::string
tempdir(const std::string& tempdir = "");
void
assert_exists(const std::string& fname);
void
assert_is_regular_file(const std::string& fname);
void
assert_dir_is_writeable(const std::string& fname);
void
assert_no_overwrite(const std::string& fname);
void
assert_can_read(const std::string& fname);
void
assert_no_cr(const std::string& infile, uint64_t nlines = 5);