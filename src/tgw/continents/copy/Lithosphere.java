package tgw.continents.copy;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import tgw.continents.copy.list.IArrayList;

import java.util.Arrays;

public class Lithosphere {

    private static final float CONTINENTAL_BASE = 1.0f;
    private static final float OCEANIC_BASE = 0.1f;
    private static final int BOOL_REGENERATE_CRUST = 1;
    public static final float SQRDMD_ROUGHNESS = 0.5f;
    public static final float SUBDUCT_RATIO = 0.5f;
    public static final float BUOYANCY_BONUS_X = 3;
    public static final int MAX_BUOYANCY_AGE = 20;
    public static final float MULINV_MAX_BUOYANCY_AGE = 1.0f / MAX_BUOYANCY_AGE;
    public static final float RESTART_ENERGY_RATIO = 0.15f;
    public static final float RESTART_SPEED_LIMIT = 2.0f;
    private static final float FLT_EPSILON = 1e-5f;
    public static final int NO_COLLISION_TIME_LIMIT = 10;
    private static final boolean DEBUG = true;
    private final int aggrOverlapAbs;
    private final float aggrOverlapRel;
    private ObjectArrayList<ObjectArrayList<PlateCollision>> collisions;
    private int cycle_count;
    private final int erosion_period;
    private final float folding_ratio;
    private final int[] globalAgemap;
    private final float[] globalHeightmap;
    /**
     * NewPlate index map of the "owner" of each map point
     */
    private int[] globalIdMap;
    private int lastCollisionCount;
    private int maxPlates;
    private final int max_cycles;
    private int numPlates;
    private double peak_Ek;
    private Plate[] plates;
    private ObjectArrayList<ObjectArrayList<PlateCollision>> subductions;
    private int tickCount;
    /**
     * Must be a power of two
     */
    private int worldSize;

    public Lithosphere(int map_side_length, float sea_level, int _erosion_period, float _folding_ratio, int aggr_ratio_abs, float aggr_ratio_rel, int num_cycles) {
        this.aggrOverlapAbs = aggr_ratio_abs;
        this.aggrOverlapRel = aggr_ratio_rel;
        this.cycle_count = 0;
        this.erosion_period = _erosion_period;
        this.folding_ratio = _folding_ratio;
        this.tickCount = 0;
        this.worldSize = map_side_length + 1;
        this.max_cycles = num_cycles;
        this.maxPlates = 0;
        this.numPlates = 0;
        final int A = this.worldSize * this.worldSize;
        float[] tmp = new float[A];
        if (Noise.sqrdmd(tmp, this.worldSize, SQRDMD_ROUGHNESS) < 0) {
            throw new IllegalArgumentException("Failed to generate height map");
        }
        float lowest = tmp[0];
        float highest = tmp[0];
        for (int i = 1; i < A; ++i) {
            lowest = Math.min(lowest, tmp[i]);
            highest = Math.max(highest, tmp[i]);
        }
        // Scale to [0 ... 1]
        for (int i = 0; i < A; ++i) {
            tmp[i] = (tmp[i] - lowest) / (highest - lowest);
        }
        float sea_threshold = 0.5f;
        float th_step = 0.5f;
        // Find the actual value in height map that produces the continent-sea
        // ratio defined be "sea_level".
        while (th_step > 0.01) {
            int count = 0;
            for (int i = 0; i < A; ++i) {
                count += tmp[i] < sea_threshold ? 1 : 0;
            }
            th_step *= 0.5f;
            if (count / (float) A < sea_level) {
                sea_threshold += th_step;
            }
            else {
                sea_threshold -= th_step;
            }
        }
        sea_level = sea_threshold;
        for (int i = 0; i < A; ++i) {
            tmp[i] = (tmp[i] > sea_level ? 1 : 0) * (tmp[i] + CONTINENTAL_BASE) + (tmp[i] <= sea_level ? 1 : 0) * OCEANIC_BASE;
        }
        // Scalp the +1 away from map side to get a power of two side length!
        // Practically only the redundant map edges become removed.
        --this.worldSize;
        this.globalHeightmap = new float[this.worldSize * this.worldSize];
        for (int i = 0; i < this.worldSize; ++i) {
            System.arraycopy(tmp, i * this.worldSize + 1, this.globalHeightmap, i * this.worldSize, this.worldSize);
        }
        this.globalIdMap = new int[this.worldSize * this.worldSize];
        this.globalAgemap = new int[this.worldSize * this.worldSize];
    }

    public static int findBound(int[] map, int length, int x0, int y0, int dx, int dy) {
        throw new AbstractMethodError();
    }

    public static int findPlate(Plate[] plates, float x, float y, int num_plates) {
        throw new AbstractMethodError();
    }

    public void createPlates(int numPlates) {
        final int mapArea = this.worldSize * this.worldSize;
        this.maxPlates = this.numPlates = numPlates;
        this.collisions = new ObjectArrayList<>(numPlates);
        this.subductions = new ObjectArrayList<>(numPlates);
        for (int i = 0; i < numPlates; ++i) {
            this.collisions.add(new ObjectArrayList<>());
            this.subductions.add(new ObjectArrayList<>());
        }
        // Initialize "Free plate center position" lookup table.
        // This way two plate centers will never be identical.
        for (int i = 0; i < mapArea; ++i) {
            this.globalIdMap[i] = i;
        }
        PlateArea[] area = Noise.fill(new PlateArea[numPlates], PlateArea::new);
        for (int i = 0; i < numPlates; ++i) {
            int plateId = this.globalIdMap[DumbMain.RANDOM.nextInt(mapArea - i)];
            int y = plateId / this.worldSize;
            int x = plateId - y * this.worldSize;
            area[i].lft = area[i].rgt = x;
            area[i].top = area[i].btm = y;
            area[i].wdt = area[i].hgt = 1;
            area[i].border.ensureCapacity(8);
            area[i].border.add(plateId); // ...and mark it as border.
            // Overwrite used entry with last unused entry in array.
            this.globalIdMap[plateId] = this.globalIdMap[mapArea - i - 1];
        }
        int[] owner = this.globalIdMap;
        Arrays.fill(owner, Integer.MAX_VALUE);
        // "Grow" plates from their origins until surface is fully populated.
        int maxBorder = 1;
        while (maxBorder != 0) {
            int i;
            for (maxBorder = i = 0; i < numPlates; ++i) {
                int bigN = area[i].border.size();
                maxBorder = Math.max(maxBorder, bigN);
                if (bigN == 0) {
                    continue;
                }
                int j = DumbMain.RANDOM.nextInt(bigN);
                int p = area[i].border.get(j);
                int cy = p / this.worldSize;
                int cx = p - cy * this.worldSize;
                int lft = cx > 0 ? cx - 1 : this.worldSize - 1;
                int rgt = cx < this.worldSize - 1 ? cx + 1 : 0;
                int top = cy > 0 ? cy - 1 : this.worldSize - 1;
                int btm = cy < this.worldSize - 1 ? cy + 1 : 0;
                int n = top * this.worldSize + cx; // North.
                int s = btm * this.worldSize + cx; // South.
                int w = cy * this.worldSize + lft; // West.
                int e = cy * this.worldSize + rgt; // East.
                if (owner[n] >= numPlates) {
                    owner[n] = i;
                    area[i].border.add(n);
                    if (area[i].top == (top + 1 & this.worldSize - 1)) {
                        area[i].top = top;
                        area[i].hgt++;
                    }
                }
                if (owner[s] >= numPlates) {
                    owner[s] = i;
                    area[i].border.add(s);
                    if (btm == (area[i].btm + 1 & this.worldSize - 1)) {
                        area[i].btm = btm;
                        area[i].hgt++;
                    }
                }
                if (owner[w] >= numPlates) {
                    owner[w] = i;
                    area[i].border.add(w);
                    if (area[i].lft == (lft + 1 & this.worldSize - 1)) {
                        area[i].lft = lft;
                        area[i].wdt++;
                    }
                }
                if (owner[e] >= numPlates) {
                    owner[e] = i;
                    area[i].border.add(e);
                    if (rgt == (area[i].rgt + 1 & this.worldSize - 1)) {
                        area[i].rgt = rgt;
                        area[i].wdt++;
                    }
                }
                // Overwrite processed point with unprocessed one.
                area[i].border.set(j, area[i].border.removeInt(area[i].border.size() - 1));
            }
        }
        this.plates = new Plate[numPlates];
        // Extract and create plates from initial terrain.
        for (int i = 0; i < numPlates; ++i) {
            area[i].wdt = area[i].wdt < this.worldSize ? area[i].wdt : this.worldSize - 1;
            area[i].hgt = area[i].hgt < this.worldSize ? area[i].hgt : this.worldSize - 1;
            final int x0 = area[i].lft;
            final int x1 = 1 + x0 + area[i].wdt;
            final int y0 = area[i].top;
            final int y1 = 1 + y0 + area[i].hgt;
            final int width = x1 - x0;
            final int height = y1 - y0;
            float[] plt = new float[width * height];
            // Copy plate's height data from global map into local map.
            for (int y = y0, j = 0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x, ++j) {
                    int k = (y & this.worldSize - 1) * this.worldSize + (x & this.worldSize - 1);
                    plt[j] = this.globalHeightmap[k] * (owner[k] == i ? 1 : 0);
                }
            }
            // Create plate.
            this.plates[i] = new Plate(plt, width, height, x0, y0, i, this.worldSize);
        }
        this.tickCount = numPlates + MAX_BUOYANCY_AGE;
        this.peak_Ek = 0;
        this.lastCollisionCount = 0;
    }

    public int[] getAgemap() {
        return this.globalAgemap;
    }

    public int getCycleCount() {
        return this.cycle_count;
    }

    public int getIterationCount() {
        return this.tickCount;
    }

    public int getPlateCount() {
        return this.numPlates;
    }

    public int[] getPlates() {
        return this.globalIdMap;
    }

    public int getSideLength() {
        return this.worldSize;
    }

    public float[] getTopography() {
        return this.globalHeightmap;
    }

    private void handlePlateCollisions(int plateId) {
        ObjectArrayList<PlateCollision> plateCollisions = this.collisions.get(plateId);
        for (int collisionId = 0, size = plateCollisions.size(); collisionId < size; ++collisionId) {
            PlateCollision collision = plateCollisions.get(collisionId);
            if (DEBUG) {
                if (plateId == collision.plateId) {
                    throw new IllegalStateException("when subducting: SRC == DEST!");
                }
            }
            Plate mainPlate = this.plates[plateId];
            Plate otherPlate = this.plates[collision.plateId];
            mainPlate.applyFriction(collision.crust);
            otherPlate.applyFriction(collision.crust);
            Plate.SegmentData mainCollisionInfo = mainPlate.getCollisionInfo(collision.x, collision.y);
            Plate.SegmentData otherCollisionInfo = otherPlate.getCollisionInfo(collision.x, collision.y);
            // Find the minimum count of collisions between two
            // continents on different plates.
            // It's minimum because large plate will get collisions
            // from all over whereas smaller plate will get just
            // a few. It's those few that matter between these two
            // plates, not what the big plate has with all the
            // other plates around it.
            int collisionCount = Math.min(mainCollisionInfo.getCollisionCount(), otherCollisionInfo.getCollisionCount());
            // Find maximum amount of collided surface area between
            // two continents on different plates.
            // Like earlier, it's the "experience" of the smaller
            // plate that matters here.
            float collisionRatio = Math.max(mainCollisionInfo.getCollisionRatio(), otherCollisionInfo.getCollisionRatio());
            if (collisionCount > this.aggrOverlapAbs || collisionRatio > this.aggrOverlapRel) {
                float amount = mainPlate.aggregateCrust(otherPlate, collision.x, collision.y);
                // Calculate new direction and speed for the
                // merged plate system, that is, for the
                // receiving plate!
                otherPlate.collide(mainPlate, collision.x, collision.y, amount);
            }
        }
        plateCollisions.clear();
    }

    private void restart() {
        final int map_area = this.worldSize * this.worldSize;
        this.cycle_count += this.max_cycles > 0 ? 1 : 0; // No increment if running for ever.
        if (this.cycle_count > this.max_cycles) {
            return;
        }
        Arrays.fill(this.globalHeightmap, 0);
        // Update height map to include all recent changes.
        for (int plateId = 0; plateId < this.numPlates; ++plateId) {
            Plate plate = this.plates[plateId];
            int x0 = (int) plate.getX();
            int y0 = (int) plate.getY();
            int x1 = x0 + plate.getWidth();
            int y1 = y0 + plate.getHeight();
            float[] heightmap = plate.getHeightmap();
            int[] agemap = plate.getAgeMap();
            // Copy first part of plate onto world map.
            for (int y = y0, index = 0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x, ++index) {
                    int wrapX = x & this.worldSize - 1;
                    int wrapY = y & this.worldSize - 1;
                    int globalIndex = wrapY * this.worldSize + wrapX;
                    float height0 = this.globalHeightmap[globalIndex];
                    float height1 = heightmap[index];
                    int age0 = this.globalAgemap[globalIndex];
                    int age1 = agemap[index];
                    this.globalAgemap[globalIndex] = (int) ((height0 * age0 + height1 * age1) / (height0 + height1));
                    this.globalHeightmap[globalIndex] += height1;
                }
            }
        }
        this.numPlates = 0;
        // create new plates IFF there are cycles left to run!
        // However, if max cycle count is "ETERNITY", then 0 < 0 + 1 always.
        if (this.cycle_count < this.max_cycles + (this.max_cycles == 0 ? 1 : 0)) {
            this.createPlates(this.numPlates = this.maxPlates);
            // Restore the ages of plates' points of crust!
            for (int plateId = 0; plateId < this.numPlates; ++plateId) {
                Plate plate = this.plates[plateId];
                final int x0 = (int) plate.getX();
                final int y0 = (int) plate.getY();
                final int x1 = x0 + plate.getWidth();
                final int y1 = y0 + plate.getHeight();
                int[] agemap = plate.getAgeMap();
                for (int y = y0, index = 0; y < y1; ++y) {
                    for (int x = x0; x < x1; ++x, ++index) {
                        int wrapX = x & this.worldSize - 1;
                        int wrapY = y & this.worldSize - 1;
                        agemap[index] = this.globalAgemap[wrapY * this.worldSize + wrapX];
                    }
                }
            }
            return;
        }
        // Add some "virginity buoyancy" to all pixels for a visual boost.
        for (int i = 0; i < map_area; ++i) {
            int crust_age = this.tickCount - this.globalAgemap[i];
            crust_age = MAX_BUOYANCY_AGE - crust_age;
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE ? 1 : 0);

            this.globalHeightmap[i] += (this.globalHeightmap[i] < CONTINENTAL_BASE ? 1 : 0) * BUOYANCY_BONUS_X * OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE;
        }
        ///////////////////////////////////////////////////////////////////////
        // This is the LAST cycle! ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        int A = (this.worldSize + 1) * (this.worldSize + 1);
        float[] tmp = new float[A];
        float[] original = new float[A];
        System.arraycopy(this.globalHeightmap, 0, original, 0, map_area);
        float h_lowest = this.globalHeightmap[0];
        float h_highest = this.globalHeightmap[0];
        for (int i = 1; i < map_area; ++i) {
            // Record elevation extremes.
            h_lowest = Math.min(h_lowest, this.globalHeightmap[i]);
            h_highest = Math.max(h_highest, this.globalHeightmap[i]);
        }
        for (int y = 0; y < this.worldSize; y += 8) {
            for (int x = 0; x < this.worldSize; x += 8) {
                int i = y * this.worldSize + x;
                this.globalHeightmap[i] = (Short.MAX_VALUE / 8) * (this.globalHeightmap[i] - h_lowest) / (h_highest - h_lowest);
                for (int j = i + 1; j < i + 8; ++j) {
                    this.globalHeightmap[j] = 0;
                }
            }
            Arrays.fill(this.globalHeightmap, (y + 1) * this.worldSize, (y + 1) * this.worldSize + this.worldSize, 0);
            Arrays.fill(this.globalHeightmap, (y + 2) * this.worldSize, (y + 2) * this.worldSize + this.worldSize, 0);
            Arrays.fill(this.globalHeightmap, (y + 3) * this.worldSize, (y + 3) * this.worldSize + this.worldSize, 0);
            for (int i = 0; i < 4; ++i) {
                this.globalHeightmap[(y + 4) * this.worldSize + i] =
                this.globalHeightmap[(y + 4) * this.worldSize + this.worldSize - i - 1] = 0;
            }
            for (int x = 4; x < this.worldSize - 4; x += 8) {
                int i = (y + 4) * this.worldSize + x;
                this.globalHeightmap[i] = (Short.MAX_VALUE / 8) * (this.globalHeightmap[i] - h_lowest) / (h_highest - h_lowest);
                for (int j = i + 1; j < i + 8; ++j) {
                    this.globalHeightmap[j] = 0;
                }
            }
            Arrays.fill(this.globalHeightmap, (y + 5) * this.worldSize, (y + 5) * this.worldSize + this.worldSize, 0);
            Arrays.fill(this.globalHeightmap, (y + 6) * this.worldSize, (y + 6) * this.worldSize + this.worldSize, 0);
            Arrays.fill(this.globalHeightmap, (y + 7) * this.worldSize, (y + 7) * this.worldSize + this.worldSize, 0);

        }
        for (int y = 0; y < this.worldSize; ++y) {
            // Copy map into fractal buffer.
            System.arraycopy(this.globalHeightmap, y * this.worldSize, tmp, y * (this.worldSize + 1), this.worldSize);
            tmp[y * (this.worldSize + 1) + this.worldSize] = this.globalHeightmap[y * this.worldSize];
        }
        // Copy last line - the one that "wraps around" the top edge.
        System.arraycopy(this.globalHeightmap, 0, tmp, this.worldSize * (this.worldSize + 1), this.worldSize);
        tmp[this.worldSize * (this.worldSize + 1) + this.worldSize] = this.globalHeightmap[0];
        // Finally create some fractal slopes!
        if (Noise.sqrdmd(tmp, this.worldSize + 1, SQRDMD_ROUGHNESS) < 0) {
            throw new IllegalArgumentException("Failed to fractalize heightmap.");
        }
        Noise.normalize(tmp, A);
        float h_range = h_highest - h_lowest;
        for (int i = 0; i < A; ++i) {
            // Restore original height range.
            tmp[i] = h_lowest + tmp[i] * h_range;
        }
        for (int y = 0; y < this.worldSize; ++y) {
            for (int x = 0; x < this.worldSize; ++x) {
                if (original[y * this.worldSize + x] > CONTINENTAL_BASE) {
                    float new_height = tmp[y * (this.worldSize + 1) + x] + original[y * this.worldSize + x] * 0.0f;
                    float alpha = (float) Math.sqrt((original[y * this.worldSize + x] - CONTINENTAL_BASE) / (h_highest - CONTINENTAL_BASE));
                    this.globalHeightmap[y * this.worldSize + x] = alpha * new_height + (1.0f - alpha) * original[y * this.worldSize + x];
                }
                else {
                    this.globalHeightmap[y * this.worldSize + x] = original[y * this.worldSize + x];
                }
            }
        }
        // Add some random noise to the map.
        Arrays.fill(tmp, 0);
        if (Noise.sqrdmd(tmp, this.worldSize + 1, SQRDMD_ROUGHNESS) < 0) {
            throw new IllegalArgumentException("Failed to generate height map again.");
        }
        Noise.normalize(tmp, A);
        // Shrink the fractal map by 1 pixel from right and bottom.
        // This makes it same size as lithosphere's height map.
        for (int i = 0; i < this.worldSize; ++i) {
            System.arraycopy(tmp, i * (this.worldSize + 1), tmp, i * this.worldSize, this.worldSize);
        }
        for (int i = 0; i < map_area; ++i) {
            if (this.globalHeightmap[i] > CONTINENTAL_BASE) {
                this.globalHeightmap[i] += tmp[i] * 2 * 0;
            }
            else {
                this.globalHeightmap[i] = 0.8f * this.globalHeightmap[i] + 0.2f * tmp[i] * CONTINENTAL_BASE;
            }
        }
        // Add a smoothing factor to sea floor.
        System.arraycopy(this.globalHeightmap, 0, original, 0, map_area);
        h_lowest = this.globalHeightmap[0];
        h_highest = this.globalHeightmap[0];
        for (int i = 1; i < map_area; ++i) {
            h_lowest = Math.min(h_lowest, this.globalHeightmap[i]);
            h_highest = Math.max(h_highest, this.globalHeightmap[i]);
        }
        for (int y = 0; y < this.worldSize; y += 4) {
            for (int x = 0; x < this.worldSize; x += 4) {
                int i = y * this.worldSize + x;
                this.globalHeightmap[i] = 4.0f * Short.MAX_VALUE * (this.globalHeightmap[i] - h_lowest) / (h_highest - h_lowest);
                for (int j = i + 1; j < i + 4; ++j) {
                    this.globalHeightmap[j] = 0;
                }
            }
            Arrays.fill(this.globalHeightmap, (y + 1) * this.worldSize, (y + 1) * this.worldSize + this.worldSize, 0);
            for (int i = 0; i < 2; ++i) {
                this.globalHeightmap[(y + 2) * this.worldSize + i] =
                this.globalHeightmap[(y + 2) * this.worldSize + this.worldSize - i - 1] = 0;
            }
            for (int x = 2; x < this.worldSize - 2; x += 4) {
                int i = (y + 2) * this.worldSize + x;
                this.globalHeightmap[i] = 4.0f * Short.MAX_VALUE * (this.globalHeightmap[i] - h_lowest) / (h_highest - h_lowest);
                for (int j = i + 1; j < i + 4; ++j) {
                    this.globalHeightmap[j] = 0;
                }
            }
            Arrays.fill(this.globalHeightmap, (y + 3) * this.worldSize, (y + 3) * this.worldSize + this.worldSize, 0);
        }

        for (int y = 0; y < this.worldSize; ++y) {
            // Copy map into fractal buffer.
            System.arraycopy(this.globalHeightmap, y * this.worldSize, tmp, y * (this.worldSize + 1), this.worldSize);
            tmp[y * (this.worldSize + 1) + this.worldSize] = this.globalHeightmap[y * this.worldSize];
        }
        // Copy last line - the one that "wraps around" the top edge.
        System.arraycopy(this.globalHeightmap, 0, tmp, this.worldSize * (this.worldSize + 1), this.worldSize);
        tmp[this.worldSize * (this.worldSize + 1) + this.worldSize] = this.globalHeightmap[0];
        // Finally create some fractal slopes!
        if (Noise.sqrdmd(tmp, this.worldSize + 1, SQRDMD_ROUGHNESS) < 0) {
            throw new IllegalArgumentException("Failed to fractalize heightmap.");
        }
        Noise.normalize(tmp, A);
        h_range = h_highest - h_lowest;
        // Restore original height range.
        for (int i = 0; i < A; ++i) {
            tmp[i] = h_lowest + tmp[i] * h_range;
        }
        for (int y = 0; y < this.worldSize; ++y) {
            for (int x = 0; x < this.worldSize; ++x) {
                if (original[y * this.worldSize + x] < CONTINENTAL_BASE) {
                    this.globalHeightmap[y * this.worldSize + x] = tmp[y * (this.worldSize + 1) + x];
                }
                else {
                    this.globalHeightmap[y * this.worldSize + x] = original[y * this.worldSize + x];
                }
            }
        }
    }

    public void update() {
        double totalVelocity = 0;
        double systemKineticEnergy = 0;
        for (int i = 0; i < this.numPlates; ++i) {
            totalVelocity += this.plates[i].getVelocity();
            systemKineticEnergy += this.plates[i].getMomentum();
        }
        if (systemKineticEnergy > this.peak_Ek) {
            this.peak_Ek = systemKineticEnergy;
        }
        // If there's no continental collisions during past iterations,
        // then interesting activity has ceased and we should restart.
        // Also if the simulation has been going on for too long already,
        // restart, because interesting stuff has most likely ended.
        if (totalVelocity < RESTART_SPEED_LIMIT || systemKineticEnergy / this.peak_Ek < RESTART_ENERGY_RATIO || this.lastCollisionCount > NO_COLLISION_TIME_LIMIT || this.tickCount > 600) {
            this.restart();
            return;
        }
        int mapArea = this.worldSize * this.worldSize;
        int[] prevIdMap = this.globalIdMap;
        this.globalIdMap = new int[mapArea];
        // Realize accumulated external forces to each plate.
        for (int plateId = 0; plateId < this.numPlates; ++plateId) {
            Plate plate = this.plates[plateId];
            plate.resetSegments();
            if (this.erosion_period > 0 && this.tickCount % this.erosion_period == 0) {
                plate.erode(CONTINENTAL_BASE);
            }
            plate.move();
        }
        int oceanic_collisions = 0;
        int continental_collisions = 0;
        // Update height and plate index maps.
        // Doing it plate by plate is much faster than doing it index wise:
        // Each plate's map's memory area is accessed sequentially and only
        // once as opposed to calculating "num_plates" indices within plate
        // maps in order to find out which plate(s) own current location.
        Arrays.fill(this.globalHeightmap, 0);
        Arrays.fill(this.globalIdMap, Integer.MAX_VALUE);
        for (int plateId = 0; plateId < this.numPlates; ++plateId) {
            Plate plate = this.plates[plateId];
            final int x0 = (int) plate.getX();
            final int y0 = (int) plate.getY();
            final int x1 = x0 + plate.getWidth();
            final int y1 = y0 + plate.getHeight();
            float[] heightmap = plate.getHeightmap();
            int[] agemap = plate.getAgeMap();
            // Copy first part of plate onto world map.
            for (int y = y0, index = 0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x, ++index) {
                    if (heightmap[index] < 2 * FLT_EPSILON) {
                        // No crust here...
                        continue;
                    }
                    int wrapX = x & this.worldSize - 1;
                    int wrapY = y & this.worldSize - 1;
                    int globalMapPos = wrapY * this.worldSize + wrapX;
                    int collidingPlateId = this.globalIdMap[globalMapPos];
                    if (collidingPlateId >= this.numPlates) {
                        // No one here yet?
                        // This plate becomes the "owner" of current location
                        // if it is the first plate to have crust on it.
                        this.globalHeightmap[globalMapPos] = heightmap[index];
                        this.globalIdMap[globalMapPos] = plateId;
                        this.globalAgemap[globalMapPos] = agemap[index];
                        continue;
                    }
                    // DO NOT ACCEPT HEIGHT EQUALITY! Equality leads to subduction
                    // of shore that's barely above sea level. It's a lot less
                    // serious problem to treat very shallow waters as continent...
                    Plate collidingPlate = this.plates[collidingPlateId];
                    int prevAge = collidingPlate.getCrustAge(wrapX, wrapY);
                    int thisAge = agemap[index];
                    // Handle subduction of oceanic crust as special case.
                    if (heightmap[index] < CONTINENTAL_BASE && (this.globalHeightmap[globalMapPos] > heightmap[index] || this.globalHeightmap[globalMapPos] + 2 * FLT_EPSILON > heightmap[index] && this.globalHeightmap[globalMapPos] < 2 * FLT_EPSILON + heightmap[index] && prevAge >= thisAge)) {
                        // This plate will be the subducting one.
                        // The level of effect that subduction has
                        // is directly related to the amount of water
                        // on top of the subducting plate.
                        float sediment = SUBDUCT_RATIO * OCEANIC_BASE * (CONTINENTAL_BASE - heightmap[index]) / CONTINENTAL_BASE;
                        // Save collision to the receiving plate's list.
                        this.subductions.get(collidingPlateId).add(new PlateCollision(plateId, wrapX, wrapY, sediment));
                        ++oceanic_collisions;
                        // Remove subducted oceanic lithosphere from plate.
                        // This is crucial for
                        // a) having correct amount of colliding crust (below)
                        // b) protecting subducted locations from receiving
                        //    crust from other subductions/collisions.
                        plate.setCrust(wrapX, wrapY, heightmap[index] - OCEANIC_BASE, thisAge);
                        if (heightmap[index] <= 0) {
                            continue; // Nothing more to collide.
                        }
                    }
                    else if (this.globalHeightmap[globalMapPos] < CONTINENTAL_BASE) {
                        float sediment = SUBDUCT_RATIO * OCEANIC_BASE * (CONTINENTAL_BASE - this.globalHeightmap[globalMapPos]) / CONTINENTAL_BASE;
                        this.subductions.get(plateId).add(new PlateCollision(collidingPlateId, wrapX, wrapY, sediment));
                        ++oceanic_collisions;
                        collidingPlate.setCrust(wrapX, wrapY, this.globalHeightmap[globalMapPos] - OCEANIC_BASE, prevAge);
                        this.globalHeightmap[globalMapPos] -= OCEANIC_BASE;
                        if (this.globalHeightmap[globalMapPos] <= 0) {
                            this.globalIdMap[globalMapPos] = plateId;
                            this.globalHeightmap[globalMapPos] = heightmap[index];
                            this.globalAgemap[globalMapPos] = agemap[index];
                            continue;
                        }
                    }
                    // Record collisions to both plates. This also creates
                    // continent segment at the collided location to plates.
                    int thisArea = plate.addCollision(wrapX, wrapY);
                    int prevArea = collidingPlate.addCollision(wrapX, wrapY);
                    // At least two plates are at same location.
                    // Move some crust from the SMALLER plate onto LARGER one.
                    if (thisArea < prevArea) {
                        float crust = heightmap[index] * this.folding_ratio;
                        // Give some...
                        this.globalHeightmap[globalMapPos] += crust;
                        collidingPlate.setCrust(wrapX, wrapY, this.globalHeightmap[globalMapPos], agemap[index]);
                        // And take some.
                        plate.setCrust(wrapX, wrapY, heightmap[index] * (1.0f - this.folding_ratio), agemap[index]);
                        // Add collision to the earlier plate's list.
                        this.collisions.get(plateId).add(new PlateCollision(collidingPlateId, wrapX, wrapY, crust));
                        ++continental_collisions;
                    }
                    else {
                        float crust = this.globalHeightmap[globalMapPos] * this.folding_ratio;
                        plate.setCrust(wrapX, wrapY, heightmap[index] + crust, this.globalAgemap[globalMapPos]);
                        collidingPlate.setCrust(wrapX, wrapY, this.globalHeightmap[globalMapPos] * (1.0f - this.folding_ratio), this.globalAgemap[globalMapPos]);
                        this.collisions.get(collidingPlateId).add(new PlateCollision(plateId, wrapX, wrapY, crust));
                        ++continental_collisions;
                        // Give the location to the larger plate.
                        this.globalHeightmap[globalMapPos] = heightmap[index];
                        this.globalIdMap[globalMapPos] = plateId;
                        this.globalAgemap[globalMapPos] = agemap[index];
                    }
                }
            }
        }
        // Update the counter of iterations since last continental collision.
        if (continental_collisions > 0) {
            ++this.lastCollisionCount;
        }
        for (int plateId = 0; plateId < this.numPlates; ++plateId) {
            ObjectArrayList<PlateCollision> plateSubductions = this.subductions.get(plateId);
            for (int j = 0, size = plateSubductions.size(); j < size; ++j) {
                PlateCollision collision = plateSubductions.get(j);
                if (DEBUG) {
                    if (plateId == collision.plateId) {
                        throw new IllegalStateException("when subducting: SRC == DEST!");
                    }
                }
                // Do not apply friction to oceanic plates.
                // This is a very cheap way to emulate slab pull.
                // Just perform subduction and on our way we go!
                this.plates[plateId].addCrustBySubduction(collision.x, collision.y, collision.crust, this.tickCount, this.plates[collision.plateId].getVelX(), this.plates[collision.plateId].getVelY());
            }
            plateSubductions.clear();
        }
        for (int plateId = 0; plateId < this.numPlates; ++plateId) {
            this.handlePlateCollisions(plateId);
        }
        int[] indexFound = new int[this.numPlates];
        // Fill divergent boundaries with new crustal material, molten magma.
        for (int y = 0, i = 0; y < BOOL_REGENERATE_CRUST * this.worldSize; ++y) {
            for (int x = 0; x < this.worldSize; ++x, ++i) {
                if (this.globalIdMap[i] >= this.numPlates) {
                    // The owner of this new crust is that neighbour plate
                    // who was located at this point before plates moved.
                    this.globalIdMap[i] = prevIdMap[i];
                    if (DEBUG) {
                        if (this.globalIdMap[i] >= this.numPlates) {
                            System.out.println("Previous index map has no owner!");
                            System.exit(1);
                        }
                    }
                    // If this is oceanic crust then add buoyancy to it.
                    // Magma that has just crystallized into oceanic crust
                    // is more buoyant than that which has had a lot of
                    // time to cool down and become denser.
                    this.globalAgemap[i] = this.tickCount;
                    this.globalHeightmap[i] = OCEANIC_BASE * BUOYANCY_BONUS_X;
                    this.plates[this.globalIdMap[i]].setCrust(x, y, OCEANIC_BASE, this.tickCount);
                }
                else if (++indexFound[this.globalIdMap[i]] != 0 && this.globalHeightmap[i] <= 0) {
                    System.out.println("Occupied point has no land mass!");
                    System.exit(1);
                }
            }
        }
        // Remove empty plates from the system.
        for (int i = 0; i < this.numPlates; ++i) {
            if (this.numPlates == 1) {
                System.out.println("ONLY ONE PLATE LEFT!");
            }
            else if (indexFound[i] == 0) {
                this.plates[i] = this.plates[this.numPlates - 1];
                indexFound[i] = indexFound[this.numPlates - 1];
                // Life is seldom as simple as seems at first.
                // Replace the moved plate's index in the index map
                // to match its current position in the array!
                for (int j = 0; j < this.worldSize * this.worldSize; ++j) {
                    if (this.globalIdMap[j] == this.numPlates - 1) {
                        this.globalIdMap[j] = i;
                    }
                }
                --this.numPlates;
                --i;
            }
        }
        // Add some "virginity buoyancy" to all pixels for a visual boost! :)
        for (int i = 0; i < mapArea; ++i) {
            // Calculate the inverted age of this piece of crust.
            // Force result to be minimum between inv. age and
            // max buoyancy bonus age.
            int crust_age = this.tickCount - this.globalAgemap[i];
            crust_age = MAX_BUOYANCY_AGE - crust_age;
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE ? 1 : 0);
            this.globalHeightmap[i] += (this.globalHeightmap[i] < CONTINENTAL_BASE ? 1 : 0) * BUOYANCY_BONUS_X * OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE;
        }
        ++this.tickCount;
    }

    public static class PlateArea {
        public IArrayList border = new IArrayList();
        public int btm;
        public int hgt;
        public int lft;
        public int rgt;
        public int top;
        public int wdt;
    }

    private static class PlateCollision {

        public float crust;
        public int plateId;
        public int x;
        public int y;

        public PlateCollision(int plateId, int x, int y, float crust) {
            this.plateId = plateId;
            this.x = x;
            this.y = y;
            this.crust = crust;
        }
    }
}
