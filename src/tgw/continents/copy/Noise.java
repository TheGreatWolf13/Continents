package tgw.continents.copy;

import java.util.function.Supplier;

public final class Noise {

    private Noise() {}

    private static float calcSum(float a, float b, float c, float d, float slope) {
        float sum = (a + b + c + d) * 0.25f;
        sum += slope * ((DumbMain.RANDOM.nextInt(Short.MAX_VALUE) << 1) - Short.MAX_VALUE);
        return sum;
    }

    public static <T> T[] fill(T[] array, Supplier<T> supplier) {
        for (int i = 0, len = array.length; i < len; i++) {
            array[i] = supplier.get();
        }
        return array;
    }

    public static void normalize(float[] arr, int size) {
        float min = arr[0];
        float max = arr[0];
        float diff;
        for (int i = 1; i < size; ++i) {
            min = Math.min(min, arr[i]);
            max = Math.max(max, arr[i]);
        }
        diff = max - min;
        if (diff > 0) {
            for (int i = 0; i < size; ++i) {
                arr[i] = (arr[i] - min) / diff;
            }
        }
    }

    private static void saveSum(int a, float[] map, float sum) {
        if ((int) map[a] == 0) {
            map[a] = sum;
        }
    }

    public static int sqrdmd(float[] map, int size, float rgh) {
        int temp = size - 1;
        /* MUST EQUAL TO 2^x + 1! */
        if ((temp & temp - 1) != 0 || (temp & 3) != 0) {
            return -1;
        }
        float slope = rgh;
        int step = size & ~1;
        /* Calculate midpoint ("diamond step"). */
        int dy = step * size;
        float sum = calcSum(map[0], map[step], map[dy], map[dy + step], slope);
        saveSum(0, map, sum);
        float center_sum = sum;
        /* Calculate each sub diamonds' center points ("square step"). */
        /* Top row. */
        int p0 = step >> 1;
        sum = calcSum(map[0], map[step], center_sum, center_sum, slope);
        saveSum(p0, map, sum);
        /* Left column. */
        int p1 = p0 * size;
        sum = calcSum(map[0], map[dy], center_sum, center_sum, slope);
        saveSum(p1, map, sum);
        final int full_size = size * size;
        map[full_size + p0 - size] = map[p0]; /* Copy top val into btm row. */
        map[p1 + size - 1] = map[p1]; /* Copy left value into right column. */
        slope *= rgh;
        step >>= 1;
        /* Enter the main loop. */
        while (step > 1) {
            //Calc midpoint of sub squares on the map ("diamond step").
            int dx = step;
            dy = step * size;
            int i = (step >> 1) * (size + 1);
            int line_jump = step * size + 1 + step - size;
            for (int y0 = 0, y1 = dy; y1 < size * size; y0 += dy, y1 += dy) {
                for (int x0 = 0, x1 = dx; x1 < size; x0 += dx, x1 += dx, i += step) {
                    sum = calcSum(map[y0 + x0], map[y0 + x1], map[y1 + x0], map[y1 + x1], slope);
                    saveSum(i, map, sum);
                }
                /* There's additional step taken at the end of last
                 * valid loop. That step actually isn't valid because
                 * the row ends right then. Thus we are forced to
                 * manually remove it after the loop so that 'i'
                 * points again to the index accessed last.
                 */
                i += line_jump - step;
            }
            /*
             * Calculate each sub diamonds' center point ("square step").
             * Diamond gets its left and right vertices from the square
             * corners of last iteration and its top and bottom vertices
             * from the "diamond step" we just performed.
             *************************************************************/
            i = step >> 1;
            p0 = step;  /* right */
            p1 = i * size + i;  /* bottom */
            /* left */
            int p2 = 0;
            /* top (wrapping edges) */
            int p3 = full_size + i - (i + 1) * size;
            /* Calculate "diamond" values for top row in map. */
            while (p0 < size) {
                sum = calcSum(map[p0], map[p1], map[p2], map[p3], slope);
                saveSum(i, map, sum);
                /* Copy it into bottom row. */
                map[full_size + i - size] = map[i];
                p0 += step;
                p1 += step;
                p2 += step;
                p3 += step;
                i += step;
            }
            /* Now that top row's values are calculated starting from
             * 'y = step >> 1' both saves us from recalculating same things
             * twice and guarantees that data will not be read beyond top
             * row of map. 'size - (step >> 1)' guarantees that data will
             * not be read beyond bottom row of map.
             */
            for (int y = step >> 1, fds = 0; y < size - (step >> 1); y += step >> 1, fds = fds != 0 ? 0 : 1) {
                p0 = step >> 1;  /* right */
                p1 = p0 * size;  /* bottom */
                p2 = -p0;  /* left */
                p3 = -p1;  /* top */
                /* For even rows add step/2. Otherwise add nothing. */
                /* Init 'x' while it's easy. */
                int x = i = p0 * fds;
                i += y * size;  /* Move 'i' into correct row. */
                p0 += i;
                p1 += i;
                /* For odd rows p2 (left) wraps around map edges. */
                p2 += i + (size - 1) * (fds != 0 ? 0 : 1);
                p3 += i;
                /* size - (step >> 1) guarantees that data will not be
                 * read beyond rightmost column of map. */
                for (; x < size - (step >> 1); x += step) {
                    sum = calcSum(map[p0], map[p1], map[p2], map[p3], slope);
                    saveSum(i, map, sum);
                    p0 += step;
                    p1 += step;
                    p2 += step;
                    p3 += step;
                    i += step;
                    /* if we start from leftmost column -> left
                     * point (p2) is going over the right border ->
                     * wrap it around into the beginning of
                     * previous rows left line. */
                    p2 -= (size - 1) * (x != 0 ? 0 : 1);
                }
                /* copy rows first element into its last */
                i = y * size;
                map[i + size - 1] = map[i];
            }
            slope *= rgh; /* reduce amount of randomness for next round */
            step >>= 1; /* split squares and diamonds in half */
        }
        return 0;
    }
}
