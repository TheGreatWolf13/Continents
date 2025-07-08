package tgw.continents;

import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongList;
import org.jetbrains.annotations.Nullable;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.random.RandomGenerator;

public final class MainTectonics {

    public static final int WORLD_SIZE = 2_048;
    private static final int NUM_PLATES = 10;
    private static final int TRENCH_WIDTH = WORLD_SIZE / 2_048;
    public static final FastRandom RANDOM = new FastRandom("Mauren".hashCode());
    public static final short[] HEIGHTMAP = new short[WORLD_SIZE * WORLD_SIZE];
    public static final int[] COLORS = new int[NUM_PLATES];
    public static final float[] BASE_CHANCE = new float[NUM_PLATES];
    public static final LongList NODES = new LongArrayList();
    private static final boolean SAVE_ALL = false;
    private static final float AREA_RATIO = 0.3f;
    public static final int CONTINENTAL_SHELF = -200;
    public static final int SUPERCONTINENT_PEAK = 1_500;
    public static final int OCEAN_BASIN = -6_400;
    private static final LongArrayList[] SUBDUCTIONS = new LongArrayList[NUM_PLATES];
    private static final LongArrayList[] COLLISIONS = new LongArrayList[NUM_PLATES];
    public static final int TRENCH_DEPTH = -12_000;
    public static final int SUBDUCTION_STEP = 1_000;
    public static final float SUBDUCTION_RATIO = 0.1f;
    public static final int MID_OCEAN_RIDGES = -2_600;
    private static final float FOLDING_RATIO = 0.001f;
    private static final float AGGR_OVERLAP_REL = 0.2F;
    //    private static final int AGGR_OVERLAP_ABS = 5_000;
    public static final boolean DEBUG = true;
    public static BufferedImage helperImage;
    public static byte[] idMap = new byte[WORLD_SIZE * WORLD_SIZE];
    private static int imageHeightmapId;
    private static int imageHelperId = 1;
    private static int imagePlateId;

    static {
        for (int i = 0; i < NUM_PLATES; i++) {
            COLORS[i] = 0xFF00_0000 | RANDOM.nextInt(0xFF_FFFF);
            BASE_CHANCE[i] = 0.2f + RANDOM.nextFloat() * 0.6f;
            //noinspection ObjectAllocationInLoop
            SUBDUCTIONS[i] = new LongArrayList();
            //noinspection ObjectAllocationInLoop
            COLLISIONS[i] = new LongArrayList();
        }
    }

    private MainTectonics() {}

    private static long createNode(int x, int y) {
        x &= 0x3fff_ffff;
        y &= 0x3fff_ffff;
        return 0b1111L << 60 | (long) x << 30 | y;
    }

    private static int fillFromNeighbour(PlateArea[] plates, RandomGenerator random, int x, int y, Direction dir) {
        int neighbourIndex = getWrappedIndex(x + dir.offX, y + dir.offY);
        byte plateId = idMap[neighbourIndex];
        if (plateId != -1) {
            return 1;
        }
        plateId = idMap[getIndex(x, y)];
        if (random.nextFloat() < BASE_CHANCE[plateId] * plates[plateId].aspectRatioChanceMult(dir)) {
            plates[plateId].extend(dir, x, y);
            idMap[neighbourIndex] = plateId;
            return 2;
        }
        return 0;
    }

    private static void fillSuperContinent() {
        Arrays.fill(HEIGHTMAP, (short) OCEAN_BASIN);
        double area = AREA_RATIO * WORLD_SIZE * WORLD_SIZE / Math.PI;
        final int radius = (int) Math.sqrt(area);
        int centerX = RANDOM.nextInt(WORLD_SIZE);
        int centerY = RANDOM.nextInt(WORLD_SIZE);
        int rSqr = radius * radius;
        for (int y = -radius; y <= radius; y++) {
            for (int x = -radius; x <= radius; x++) {
                int dr = x * x + y * y;
                if (dr <= rSqr) {
                    HEIGHTMAP[getWrappedIndex(centerX + x, centerY + y)] = (short) (CONTINENTAL_SHELF + (SUPERCONTINENT_PEAK - CONTINENTAL_SHELF) * (1.0 - Math.sqrt(dr) / radius));
                }
            }
        }
    }

    public static int getHeightmapColor(short height) {
        if (height == Short.MIN_VALUE) {
            return 0xFFFF_00DC;
        }
        if (height >= 0) {
            if (height < 100) {
                return 0xFF00_7575;
            }
            if (height < 200) {
                return 0xFF34_9084;
            }
            if (height < 500) {
                return 0xFF6C_AB93;
            }
            if (height < 1_000) {
                return 0xFFA9_CAA4;
            }
            if (height < 1_500) {
                return 0xFFE5_DDB4;
            }
            if (height < 2_000) {
                return 0xFFC2_977C;
            }
            if (height < 3_000) {
                return 0xFF9E_5447;
            }
            if (height < 4_000) {
                return 0xFFAC_ACAC;
            }
            return 0xFFF2_EFEE;
        }
        if (height >= -200) {
            return 0xFFAB_E1FD;
        }
        if (height >= -2_000) {
            return 0xFF91_C9EF;
        }
        if (height >= -3_000) {
            return 0xFF77_B6DC;
        }
        if (height >= -4_000) {
            return 0xFF5F_A1CA;
        }
        if (height >= -5_000) {
            return 0xFF45_8FBB;
        }
        if (height >= -6_000) {
            return 0xFF30_7BAC;
        }
        if (height >= -7_000) {
            return 0xFF1A_689A;
        }
        return 0xFF13_4E70;
    }

    private static int getIndex(int x, int y) {
        return x + y * WORLD_SIZE;
    }

    private static int getNodeX(long node) {
        return (int) (node >> 30 & 0x3fff_ffff);
    }

    private static int getNodeY(long node) {
        return (int) (node & 0x3fff_ffff);
    }

    private static Direction getRandomDir(long node, RandomGenerator random) {
        int size = Long.bitCount(node & 0b1111L << 60);
        int dirIndex = random.nextInt(size);
        int count = 0;
        for (int i = 0; ; ++i) {
            if ((node & 1L << i + 60) != 0) {
                if (count++ == dirIndex) {
                    return Direction.VALUES[i];
                }
            }
        }
    }

    private static int getWrappedIndex(int x, int y) {
        return getIndex(wrap(x), wrap(y));
    }

    private static void handleCollisions(NewPlate[] plates, int plateId) {
        NewPlate plate = plates[plateId];
        LongArrayList collisions = COLLISIONS[plateId];
        for (int i = 0, len = collisions.size(); i < len; ++i) {
            long collision = collisions.getLong(i);
            NewPlate collidingPlate = plates[unpackCollisionPlate(collision)];
            short crust = unpackCollisionCrust(collision);
            plate.applyFriction(crust);
            collidingPlate.applyFriction(crust);
            int collisionX = unpackCollisionX(collision);
            int collisionY = unpackCollisionY(collision);
            NewPlate.SegmentData mainCollisionInfo = plate.getCollisionInfo(collisionX, collisionY);
            NewPlate.SegmentData otherCollisionInfo = collidingPlate.getCollisionInfo(collisionX, collisionY);
            int collisionCount = Math.min(mainCollisionInfo.getCollisionCount(), otherCollisionInfo.getCollisionCount());
            float collisionRatio = Math.max(mainCollisionInfo.getCollisionRatio(), otherCollisionInfo.getCollisionRatio());
            if (/*collisionCount > AGGR_OVERLAP_ABS ||*/ collisionRatio > AGGR_OVERLAP_REL) {
                short aggregatedCrust = plate.aggregateCrust(collidingPlate, collisionX, collisionY);
                collidingPlate.collide(plate, collisionX, collisionY, aggregatedCrust);
            }
        }
        collisions.clear();
    }

    private static void handleSubductions(NewPlate[] plates, int plateId) {
        NewPlate plate = plates[plateId];
        LongArrayList subductions = SUBDUCTIONS[plateId];
        for (int i = 0, len = subductions.size(); i < len; ++i) {
            long collision = subductions.getLong(i);
            NewPlate collidingPlate = plates[unpackCollisionPlate(collision)];
            int globalX = unpackCollisionX(collision);
            int globalY = unpackCollisionY(collision);
            plate.addCrustBySubduction(globalX, globalY, unpackCollisionCrust(collision), collidingPlate.getVelX(), collidingPlate.getVelY());
            collidingPlate.makeSubdutionRamp(globalX, globalY);
        }
        subductions.clear();
    }

    public static void main(String[] args) throws IOException {
        helperImage = new BufferedImage(WORLD_SIZE, WORLD_SIZE, BufferedImage.TYPE_INT_ARGB);
        Arrays.fill(idMap, (byte) -1);
        PlateArea[] plateAreas = new PlateArea[NUM_PLATES];
        for (int plateId = 0; plateId < NUM_PLATES; plateId++) {
            byte oldId = 0;
            int x = -1;
            int y = -1;
            int index = -1;
            while (oldId != -1) {
                x = RANDOM.nextInt(WORLD_SIZE);
                y = RANDOM.nextInt(WORLD_SIZE);
                index = getIndex(x, y);
                oldId = idMap[index];
            }
            idMap[index] = (byte) plateId;
            NODES.add(createNode(x, y));
            //noinspection ObjectAllocationInLoop
            plateAreas[plateId] = new PlateArea(x, y);
        }
        int count = 0;
        while (!NODES.isEmpty()) {
            int size = NODES.size();
            for (int i = 0; i < size; i++) {
                long node = NODES.getLong(i);
                Direction dir = getRandomDir(node, RANDOM);
                int nodeX = getNodeX(node);
                int nodeY = getNodeY(node);
                int result = fillFromNeighbour(plateAreas, RANDOM, nodeX, nodeY, dir);
                if (result > 0) {
                    node = removeDir(node, dir);
                    if (result > 1) {
                        NODES.add(removeDir(createNode(wrap(nodeX + dir.offX), wrap(nodeY + dir.offY)), dir.opposite()));
                    }
                    if (shouldNodeBeRemoved(node)) {
                        NODES.remove(i--);
                        --size;
                    }
                    else {
                        NODES.set(i, node);
                    }
                }
            }
            if ((count++ & 15) == 0) {
                savePlates(null, false);
            }
        }
        fillSuperContinent();
        NewPlate[] plates = new NewPlate[NUM_PLATES];
        for (int i = 0; i < NUM_PLATES; ++i) {
            PlateArea area = plateAreas[i];
            //noinspection ObjectAllocationInLoop
            plates[i] = new NewPlate(HEIGHTMAP, area.x0, area.y0, area.width, area.height, idMap, (byte) i);
        }
        Arrays.sort(plates, Comparator.comparingDouble(NewPlate::getMass));
        Arrays.fill(idMap, (byte) -1);
        for (int plateId = 0; plateId < NUM_PLATES; plateId++) {
            NewPlate plate = plates[plateId];
            short[] heightmap = plate.getHeightmap();
            int x0 = (int) plate.getX();
            int y0 = (int) plate.getY();
            int x1 = x0 + plate.getWidth();
            int y1 = y0 + plate.getHeight();
            for (int y = y0, index = 0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x, ++index) {
                    if (heightmap[index] == Short.MIN_VALUE) {
                        //No crust here
                        continue;
                    }
                    int globalX = x & WORLD_SIZE - 1;
                    int globalY = y & WORLD_SIZE - 1;
                    int globalIndex = globalY * WORLD_SIZE + globalX;
                    idMap[globalIndex] = (byte) plateId;
                }
            }
            plate.save(plateId);
        }
        saveHeightmap();
        savePlates(null, true);
//        for (int plateId = 0; plateId < NUM_PLATES; plateId++) {
//            NewPlate plate = plates[plateId];
//            System.out.printf("Extending plate %3d\n", plateId);
//            plate.extendPlate(0, 0);
//            plate.save(plateId);
//        }
        for (int tick = 0; tick < 10; tick++) {
            for (int i = 0; i < WORLD_SIZE; ++i) {
                for (int j = 0; j < WORLD_SIZE; ++j) {
                    helperImage.setRGB(i, j, 0);
                }
            }
            mainLoop(plates);
            saveHeightmap();
            savePlates(plates, true);
            ImageIO.write(helperImage, "png", new File(String.format(Location.ROOT_FOLDER + "\\helper_%03d.png", imageHelperId++)));
        }
    }

    private static void mainLoop(NewPlate[] plates) {
        //Erode and move
        for (int plateId = 0; plateId < NUM_PLATES; ++plateId) {
            NewPlate plate = plates[plateId];
            plate.resetSegments();
            //TODO Erosion
            plate.move();
        }
        //Update global maps, recording collisions
        Arrays.fill(HEIGHTMAP, Short.MIN_VALUE);
        Arrays.fill(idMap, (byte) -1);
        int missingCrust = WORLD_SIZE * WORLD_SIZE;
        for (int plateId = 0; plateId < NUM_PLATES; ++plateId) {
            NewPlate plate = plates[plateId];
            int x0 = (int) plate.getX();
            int y0 = (int) plate.getY();
            int x1 = x0 + plate.getWidth();
            int y1 = y0 + plate.getHeight();
            short[] heightmap = plate.getHeightmap();
            for (int y = y0, index = 0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x, ++index) {
                    if (heightmap[index] == Short.MIN_VALUE) {
                        //No crust here
                        continue;
                    }
                    int globalX = x & WORLD_SIZE - 1;
                    int globalY = y & WORLD_SIZE - 1;
                    int globalIndex = globalY * WORLD_SIZE + globalX;
                    int collidingPlateId = idMap[globalIndex];
                    if (collidingPlateId == -1) {
                        //Empty space
                        HEIGHTMAP[globalIndex] = heightmap[index];
                        idMap[globalIndex] = (byte) plateId;
                        --missingCrust;
                        continue;
                    }
                    short currentHeight = HEIGHTMAP[globalIndex];
                    short evaluatingHeight = heightmap[index];
                    NewPlate collidingPlate = plates[collidingPlateId];
                    if (evaluatingHeight < CONTINENTAL_SHELF && currentHeight >= evaluatingHeight) {
                        //Subduct this plate
                        boolean fullySubmerged = true;
                        testForFullySubmerged:
                        for (int dx = -1; dx <= 1; dx++) {
                            int newX = globalX + dx & WORLD_SIZE - 1;
                            for (int dy = -1; dy <= 1; dy++) {
                                int newY = globalY + dy & WORLD_SIZE - 1;
                                int newIndex = newY * WORLD_SIZE + newX;
                                if (idMap[newIndex] != collidingPlateId) {
                                    fullySubmerged = false;
                                    break testForFullySubmerged;
                                }
                            }
                        }
                        if (fullySubmerged) {
                            short sediment = (short) (SUBDUCTION_RATIO * (CONTINENTAL_SHELF - evaluatingHeight));
                            SUBDUCTIONS[collidingPlateId].add(packCollision(plateId, globalX, globalY, sediment));
                            plate.setCrust(globalX, globalY, (short) (evaluatingHeight - SUBDUCTION_STEP));
                            helperImage.setRGB(globalX, globalY, 0xFF00_00FF);
                        }
                        continue;
                    }
                    if (currentHeight < CONTINENTAL_SHELF) {
                        //Subduct colliding plate
                        short sediment = (short) (SUBDUCTION_RATIO * (CONTINENTAL_SHELF - currentHeight));
                        SUBDUCTIONS[plateId].add(packCollision(collidingPlateId, globalX, globalY, sediment));
                        collidingPlate.setCrust(globalX, globalY, (short) (currentHeight - SUBDUCTION_STEP));
                        HEIGHTMAP[globalIndex] -= SUBDUCTION_STEP;
                        if (HEIGHTMAP[globalIndex] < TRENCH_DEPTH) {
                            idMap[globalIndex] = (byte) plateId;
                            HEIGHTMAP[globalIndex] = evaluatingHeight;
                        }
                        helperImage.setRGB(globalX, globalY, 0xFF00_00FF);
                        continue;
                    }
                    //Continental collisions
                    int thisArea = plate.addCollision(globalX, globalY);
                    int otherArea = collidingPlate.addCollision(globalX, globalY);
                    if (thisArea < otherArea) {
                        short crust = (short) ((evaluatingHeight - CONTINENTAL_SHELF) * FOLDING_RATIO);
                        HEIGHTMAP[globalIndex] += crust;
                        collidingPlate.setCrust(globalX, globalY, (short) (currentHeight + crust));
                        plate.setCrust(globalX, globalY, (short) (evaluatingHeight - crust));
                        COLLISIONS[plateId].add(packCollision(collidingPlateId, globalX, globalY, crust));
                    }
                    else {
                        short crust = (short) ((currentHeight - CONTINENTAL_SHELF) * FOLDING_RATIO);
                        plate.setCrust(globalX, globalY, (short) (evaluatingHeight + crust));
                        collidingPlate.setCrust(globalX, globalY, (short) (currentHeight - crust));
                        COLLISIONS[collidingPlateId].add(packCollision(plateId, globalX, globalY, crust));
                        HEIGHTMAP[globalIndex] = heightmap[index];
                        idMap[globalIndex] = (byte) plateId;
                    }
                }
            }
        }
        //Handle Subductions
        for (int plateId = 0; plateId < NUM_PLATES; ++plateId) {
            handleSubductions(plates, plateId);
        }
        //Handle Collisions
        for (int plateId = 0; plateId < NUM_PLATES; ++plateId) {
            //This affects motion and "glues" plates together
            handleCollisions(plates, plateId);
        }
        //Reform crust
        while (missingCrust > 0) {
            for (int y = 0, index = 0; y < WORLD_SIZE; ++y) {
                for (int x = 0; x < WORLD_SIZE; ++x, ++index) {
                    if (idMap[index] == -1) {
                        for (Direction dir : Direction.getRandomOrder(RANDOM)) {
                            int sideX = x + dir.offX & WORLD_SIZE - 1;
                            int sideY = y + dir.offY & WORLD_SIZE - 1;
                            int sideIndex = sideY * WORLD_SIZE + sideX;
                            byte sidePlateId = idMap[sideIndex];
                            if (sidePlateId != -1) {
                                --missingCrust;
                                idMap[index] = sidePlateId;
                                HEIGHTMAP[index] = MID_OCEAN_RIDGES;
                                plates[sidePlateId].setCrust(x, y, (short) MID_OCEAN_RIDGES);
                                helperImage.setRGB(x, y, 0xFFFF_FFFF);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    private static long packCollision(int plateId, int x, int y, short crust) {
        return (plateId & 0xFFL) << 48 | (crust & 0xFFFFL) << 32 | (x & 0xFFFFL) << 16 | y & 0xFFFFL;
    }

    private static long removeDir(long node, Direction dir) {
        long mask = 1L << 60 + dir.ordinal();
        return node & ~mask;
    }

    private static void saveHeightmap() throws IOException {
        BufferedImage image = new BufferedImage(WORLD_SIZE, WORLD_SIZE, BufferedImage.TYPE_INT_ARGB);
        for (int y = 0; y < WORLD_SIZE; ++y) {
            for (int x = 0; x < WORLD_SIZE; ++x) {
                image.setRGB(x, y, getHeightmapColor(HEIGHTMAP[getIndex(x, y)]));
            }
        }
        ImageIO.write(image, "png", new File(String.format(Location.ROOT_FOLDER + "\\heightmap_%03d.png", imageHeightmapId++)));
    }

    private static void savePlates(NewPlate @Nullable [] plates, boolean force) throws IOException {
        if (!SAVE_ALL && !force) {
            return;
        }
        BufferedImage image = new BufferedImage(WORLD_SIZE, WORLD_SIZE, BufferedImage.TYPE_INT_ARGB);
        for (int y = 0; y < WORLD_SIZE; ++y) {
            for (int x = 0; x < WORLD_SIZE; ++x) {
                byte plateId = idMap[getIndex(x, y)];
                image.setRGB(x, y, plateId == -1 ? 0 : COLORS[plateId]);
            }
        }
        if (plates != null) {
            for (NewPlate plate : plates) {
                plate.attachDebugInfo(image);
            }
        }
        ImageIO.write(image, "png", new File(String.format(Location.ROOT_FOLDER + "\\plates_%03d.png", imagePlateId++)));
    }

    private static boolean shouldNodeBeRemoved(long node) {
        return (node & 0b1111L << 60) == 0;
    }

    private static short unpackCollisionCrust(long collision) {
        return (short) (collision >> 32 & 0xFFFF);
    }

    private static byte unpackCollisionPlate(long collision) {
        return (byte) (collision >> 48 & 0xFF);
    }

    private static int unpackCollisionX(long collision) {
        return (int) (collision >> 16 & 0xFFFF);
    }

    private static int unpackCollisionY(long collision) {
        return (int) (collision & 0xFFFF);
    }

    public static int wrap(int value) {
        return value & WORLD_SIZE - 1;
    }

    private static class PlateArea {

        int height;
        final int originX;
        final int originY;
        int width;
        int x0;
        int x1;
        int y0;
        int y1;

        public PlateArea(int originX, int originY) {
            this.originX = originX;
            this.originY = originY;
            this.x0 = originX;
            this.x1 = originX;
            this.y0 = originY;
            this.y1 = originY;
            this.width = 1;
            this.height = 1;
        }

        public float aspectRatioChanceMult(Direction dir) {
            return switch (dir) {
                case SOUTH, NORTH -> {
                    float aspectRatio = (float) this.width / this.height;
                    yield aspectRatio * aspectRatio;
                }
                case EAST, WEST -> {
                    float aspectRatio = (float) this.height / this.width;
                    yield aspectRatio * aspectRatio;
                }
            };
        }

        public void extend(Direction dir, int x, int y) {
            switch (dir) {
                case NORTH -> {
                    if (this.y0 == y) {
                        this.y0 = wrap(this.y0 - 1);
                        if (++this.height > WORLD_SIZE) {
                            this.height = WORLD_SIZE;
                        }
                    }
                }
                case SOUTH -> {
                    if (this.y1 == y) {
                        this.y1 = wrap(this.y1 + 1);
                        if (++this.height > WORLD_SIZE) {
                            this.height = WORLD_SIZE;
                        }
                    }

                }
                case EAST -> {
                    if (this.x1 == x) {
                        this.x1 = wrap(this.x1 + 1);
                        if (++this.width > WORLD_SIZE) {
                            this.width = WORLD_SIZE;
                        }
                    }
                }
                case WEST -> {
                    if (this.x0 == x) {
                        this.x0 = wrap(this.x0 - 1);
                        if (++this.width > WORLD_SIZE) {
                            this.width = WORLD_SIZE;
                        }
                    }
                }
            }
        }
    }

    public enum Direction {
        NORTH(0, -1, 1),
        SOUTH(0, 1, 0),
        EAST(1, 0, 3),
        WEST(-1, 0, 2);

        public static final Direction[] VALUES = values();
        private static final Direction[] RANDOM = values();
        public final int offX;
        public final int offY;
        private final int oppositeIndex;

        Direction(int offX, int offY, int oppositeIndex) {
            this.offX = offX;
            this.offY = offY;
            this.oppositeIndex = oppositeIndex;
        }

        public static Direction[] getRandomOrder(RandomGenerator random) {
            int index = random.nextInt(4);
            Direction temp = RANDOM[index];
            RANDOM[index] = RANDOM[3];
            RANDOM[3] = temp;
            index = random.nextInt(3);
            temp = RANDOM[index];
            RANDOM[index] = RANDOM[2];
            RANDOM[2] = temp;
            index = random.nextInt(2);
            temp = RANDOM[index];
            RANDOM[index] = RANDOM[1];
            RANDOM[1] = temp;
            return RANDOM;
        }

        public Direction opposite() {
            return VALUES[this.oppositeIndex];
        }
    }
}
