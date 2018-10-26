/*************************************************************************
zdvis: Lagrangian Visualization for Vector, Tensor, and Multifield Data.

Author: Zi'ang Ding

Copyright (c) 2016-2018, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/
#ifndef _ZD_LIB_FILE_TOOL_HPP_
#define _ZD_LIB_FILE_TOOL_HPP_

#include <string>
#include "../define.hpp"
#include "../configure.hpp"

namespace ZD {
    class CFileTool {
    public:
        /* get path from pathname */
        static std::string GetFilePath(std::string &file)
        {
            return GetFilePath(file.c_str());
        }
        static std::string GetFilePath(const char *file)
        {
            std::string str = file;
            size_t pos = str.find_last_of(SLASH);
            return str.substr(0, pos + 1);
        }

        /* get name from pathname */
        static std::string GetFileName(std::string& file)
        {
            return GetFileName(file.c_str());
        }
        static std::string GetFileName(const char *file)
        {
            std::string str = file;
            size_t pos_left = str.find_last_of(SLASH);
            size_t pos_right = str.find_last_of('.');
            return str.substr(pos_left+1, pos_right-pos_left-1);
        }

        /* get file extension from pathname */
        static std::string GetFileExtension(std::string &file)
        {
            return GetFileExtension(file.c_str());
        }
        static std::string GetFileExtension(const char *file)
        {
            std::string str = file;
            size_t pos = str.find_last_of('.');
            return str.substr(pos + 1, str.size() - pos - 1);
        }

        /* remove extension */
        static std::string RemoveFileExtension(std::string &file)
        {
            return RemoveFileExtension(file.c_str());
        }
        static std::string RemoveFileExtension(const char *file)
        {
            std::string str = file;
            size_t pos = str.find_last_of('.');
            return str.substr(0, pos);
        }

        /* read line function */
        static void ReadOneLine(FILE *fp, char *line)
        {
            int id = 0;
            char c = fgetc(fp);
            while (c != '\n' && c != EOF) {
                line[id] = c;
                id++;
                c = fgetc(fp);
            }
            line[id] = '\0';
        }

        /* skip lines function */
        static void SkipLines(FILE *fp, const int k) 
        {
            for (int i = 0; i < k; ++i) {
                char buffer[1024];
                ReadOneLine(fp, buffer);
            }
        }
    };
}

#endif
