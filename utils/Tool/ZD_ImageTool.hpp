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
#ifndef _ZD_LIB_IMAGE_TOOL_HPP_
#define _ZD_LIB_IMAGE_TOOL_HPP_

#include <stdio.h>

#include "../Base/ZD_Point.hpp"

// TEEM
#include <teem/nrrd.h>
#include <teem/biff.h>

#include "../define.hpp"


namespace ZD {
    class CImageTool {
    public:
        static int SaveAsPNG(unsigned char *pixels, const int width,
            const int height, const char *pathname, FILE *fp = stdout)
        {
            if (pathname == NULL || pixels == NULL || width == 0 || height == 0)
                return -1;

            Nrrd *dst = nrrdNew();
            size_t ss[3];
            ss[0] = 3;
            ss[1] = width;
            ss[2] = height;

            nrrdAlloc_nva(dst, nrrdTypeUChar, 3, ss);
            unsigned char *data = (unsigned char *)(dst->data);
            int size = 3 * width * height;
            memcpy(data, pixels, sizeof(unsigned char)*size);

            nrrdAxisInfoSet_va(dst, nrrdAxisInfoSpacing, 1.0, 1.0, 1.0);
            nrrdAxisInfoSet_va(dst, nrrdAxisInfoKind, nrrdKindSpace, nrrdKindSpace, nrrdKindSpace);
            nrrdAxisInfoSet_va(dst, nrrdAxisInfoUnits, "mm", "mm", "mm");

            if (nrrdSave(pathname, dst, NULL)) {
                fprintf(fp, "cannot save file %s:\n%s\n", pathname, biffGetDone(NRRD));
                return -1;
            }

            nrrdNuke(dst);

            return 0;
        }

        static int SaveAsPNG(CPoint<double, 1> *pixels, const int width,
            const int height, const char *pathname, FILE *fp = stdout)
        {
            if (pathname == NULL || pixels == NULL || width == 0 || height == 0)
                return -1;

            float min_v = pixels[0][0];
            float max_v = pixels[0][0];
            for (int i = 1; i < width*height; ++i) {
                if (pixels[i][0] > max_v)
                    max_v = pixels[i][0];
                if (pixels[i][0] < min_v)
                    min_v = pixels[i][0];
            }
            unsigned char *rgb = new unsigned char[width*height * 3];
            for (int i = 0; i < width*height; ++i) {
                rgb[i * 3 + 0] = (unsigned char)((pixels[i][0] - min_v) / (max_v - min_v) * 255.0);
                rgb[i * 3 + 2] = rgb[i * 3 + 1] = rgb[i * 3 + 0];
            }

            SaveAsPNG(rgb, width, height, pathname, fp);

            SafeDeleteArray(rgb);

            return 0;
        }

        static int SaveAsPNG(CPoint<double, 3> *pixels, const int width,
            const int height, const char *pathname, FILE *fp = stdout)
        {
            if (pathname == NULL || pixels == NULL || width == 0 || height == 0)
                return -1;

            unsigned char *rgb = new unsigned char[width*height * 3];
            for (int i = 0; i < width*height; ++i) {
                rgb[i * 3 + 0] = (unsigned char)(pixels[i][0] * 255.0);
                rgb[i * 3 + 1] = (unsigned char)(pixels[i][1] * 255.0);
                rgb[i * 3 + 2] = (unsigned char)(pixels[i][2] * 255.0);
            }

            SaveAsPNG(rgb, width, height, pathname, fp);

            SafeDeleteArray(rgb);

            return 0;
        }

        static int SaveAsPNG(double *pixels, const int width,
            const int height, const char *pathname, FILE *fp = stdout)
        {
            if (pathname == NULL || pixels == NULL || width == 0 || height == 0)
                return -1;

            double min_v = pixels[0];
            double max_v = pixels[0];
            for (int i = 1; i < width*height; ++i) {
                if (pixels[i] > max_v)
                    max_v = pixels[i];
                if (pixels[i] < min_v)
                    min_v = pixels[i];
            }
            unsigned char *rgb = new unsigned char[width*height * 3];
            for (int i = 0; i < width*height; ++i) {
                rgb[i * 3 + 0] = (unsigned char)((pixels[i] - min_v) / (max_v - min_v) * 255.0);
                rgb[i * 3 + 2] = rgb[i * 3 + 1] = rgb[i * 3 + 0];
            }

            SaveAsPNG(rgb, width, height, pathname, fp);

            SafeDeleteArray(rgb);

            return 0;
        }
    };
}

#endif
