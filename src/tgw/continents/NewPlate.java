package tgw.continents;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import tgw.continents.copy.DumbMain;
import tgw.continents.copy.Noise;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class NewPlate {
    private static final float DEFORMATION_COEF = 2.0f;
    private double accX;
    private double accY;
    private int activeContinent;
    private double centerOfMassX;
    private double centerOfMassY;
    private int height;
    private short[] heightmap;
    private int imageId;
    private double mass;
    private final int rotationDir;
    private final ObjectList<SegmentData> segmentData;
    private int[] segmentIdMap;
    private double speed;
    private double velX;
    private double velY;
    private int width;
    private double x;
    private double y;

    public NewPlate(short[] heightmap, int x, int y, int width, int height, byte[] idMap, byte id) {
        this.segmentData = new ObjectArrayList<>();
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;
        int area = this.width * this.height;
        this.heightmap = new short[area];
        this.segmentIdMap = new int[area];
        Arrays.fill(this.heightmap, Short.MIN_VALUE);
        Arrays.fill(this.segmentIdMap, Integer.MAX_VALUE);
        this.speed = MainTectonics.RANDOM.nextDouble() * 1.5;
        this.rotationDir = MainTectonics.RANDOM.nextBoolean() ? 1 : -1;
        double angle = 2 * Math.PI * MainTectonics.RANDOM.nextDouble();
        this.velX = Math.cos(angle);
        this.velY = Math.sin(angle);
        int y1 = y + height;
        int x1 = x + width;
        for (int dy = y, index = 0; dy < y1; ++dy) {
            for (int dx = x; dx < x1; ++dx, ++index) {
                int globalIndex = (dy & MainTectonics.WORLD_SIZE - 1) * MainTectonics.WORLD_SIZE + (dx & MainTectonics.WORLD_SIZE - 1);
                if (idMap[globalIndex] == id) {
                    short crust = heightmap[globalIndex];
                    this.heightmap[index] = crust;
                    int localMass = crust >= 0 ? crust : -2 * crust;
                    this.mass += localMass;
                    this.centerOfMassX += localMass * (dx - x);
                    this.centerOfMassY += localMass * (dy - y);
                }
            }
        }
        this.centerOfMassX /= this.mass;
        this.centerOfMassY /= this.mass;
    }

    public int addCollision(int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        int segmentId = this.segmentIdMap[index];
        if (segmentId >= this.segmentData.size()) {
            segmentId = this.createSegment(localX, localY);
        }
        SegmentData segmentData = this.segmentData.get(segmentId);
        ++segmentData.collisionCount;
        return segmentData.area;
    }

    private void addCrustByCollision(int globalX, int globalY, short crust) {
        short crustAmount = this.getCrustAmount(globalX, globalY);
        if (crustAmount == Short.MIN_VALUE) {
            crustAmount = 0;
        }
        this.setCrust(globalX, globalY, (short) (crustAmount + crust));
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        this.segmentIdMap[index] = this.activeContinent;
        SegmentData data = this.segmentData.get(this.activeContinent);
        ++data.area;
        if (localY < data.y0) {
            data.y0 = localY;
        }
        if (localY > data.y1) {
            data.y1 = localY;
        }
        if (localX < data.x0) {
            data.x0 = localX;
        }
        if (localX > data.x1) {
            data.x1 = localX;
        }
    }

    public void addCrustBySubduction(int globalX, int globalY, short crust, double subVelX, double subVelY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        if (this.velX * subVelX + this.velY * subVelY > 0) {
            subVelX -= this.velX;
            subVelY -= this.velY;
        }
        float offset = DumbMain.RANDOM.nextFloat();
        offset *= offset * offset * (2 * DumbMain.RANDOM.nextInt(2) - 1);
        subVelX = 10 * subVelX + 3 * offset;
        subVelY = 10 * subVelY + 3 * offset;
        localX += (int) subVelX;
        localY += (int) subVelY;
        if (this.width == MainTectonics.WORLD_SIZE) {
            localX &= this.width - 1;
        }
        if (this.height == MainTectonics.WORLD_SIZE) {
            localY &= this.height - 1;
        }
        index = localY * this.width + localX;
        if (index >= 0 && index < this.width * this.height && this.heightmap[index] != Short.MIN_VALUE) {
            MainTectonics.helperImage.setRGB(localX + (int) this.x & MainTectonics.WORLD_SIZE - 1, localY + (int) this.y & MainTectonics.WORLD_SIZE - 1, 0xFF00_77FF);
            this.heightmap[index] += crust;
            this.mass += crust;
        }
    }

    public short aggregateCrust(NewPlate destinationPlate, int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        int segmentId = this.segmentIdMap[index];
        SegmentData segmentData = this.segmentData.get(segmentId);
        if (segmentData.area == 0) {
            return 0;
        }
        destinationPlate.selectCollisionSegment(globalX, globalY);
        // Wrap coordinates around world edges to safeguard subtractions.
        globalX += MainTectonics.WORLD_SIZE;
        globalY += MainTectonics.WORLD_SIZE;
        double oldMass = this.mass;
        // Add all the collided continent's crust to destination plate.
        for (int y = segmentData.y0; y <= segmentData.y1; ++y) {
            for (int x = segmentData.x0; x <= segmentData.x1; ++x) {
                int i = y * this.width + x;
                if (this.segmentIdMap[i] == segmentId && this.heightmap[i] > MainTectonics.CONTINENTAL_SHELF) {
                    destinationPlate.addCrustByCollision(globalX + x - localX & MainTectonics.WORLD_SIZE - 1, globalY + y - localY & MainTectonics.WORLD_SIZE - 1, this.heightmap[i]);
                    this.mass -= this.heightmap[i];
                    this.heightmap[i] = Short.MIN_VALUE;
                }
            }
        }
        segmentData.area = 0; // Mark segment as non-existent.
        return (short) (oldMass - this.mass);
    }

    public void applyFriction(short crust) {
        if (this.mass > 0) {
            double velDec = DEFORMATION_COEF * crust / this.mass;
            this.speed -= Math.min(velDec, this.speed);
        }
    }

    public void attachDebugInfo(BufferedImage image) {
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                image.setRGB((int) (this.centerOfMassX + this.x + dx) & MainTectonics.WORLD_SIZE - 1, (int) (this.centerOfMassY + this.y + dy) & MainTectonics.WORLD_SIZE - 1, 0xFFFF_0000);
            }
        }
        int vx = (int) (this.speed * this.velX * 50);
        int vy = (int) (this.speed * this.velY * 50);
        if (vx >= 0) {
            for (int dx = 0; dx <= vx; ++dx) {
                image.setRGB((int) (this.centerOfMassX + this.x + dx) & MainTectonics.WORLD_SIZE - 1, (int) (this.centerOfMassY + this.y) & MainTectonics.WORLD_SIZE - 1, 0xFF00_0000);
            }
        }
        else {
            for (int dx = vx; dx <= 0; ++dx) {
                image.setRGB((int) (this.centerOfMassX + this.x + dx) & MainTectonics.WORLD_SIZE - 1, (int) (this.centerOfMassY + this.y) & MainTectonics.WORLD_SIZE - 1, 0xFF00_0000);
            }
        }
        if (vy >= 0) {
            for (int dy = 0; dy <= vy; ++dy) {
                image.setRGB((int) (this.centerOfMassX + this.x) & MainTectonics.WORLD_SIZE - 1, (int) (this.centerOfMassY + this.y + dy) & MainTectonics.WORLD_SIZE - 1, 0xFF00_0000);
            }
        }
        else {
            for (int dy = vy; dy <= 0; ++dy) {
                image.setRGB((int) (this.centerOfMassX + this.x) & MainTectonics.WORLD_SIZE - 1, (int) (this.centerOfMassY + this.y + dy) & MainTectonics.WORLD_SIZE - 1, 0xFF00_0000);
            }
        }
    }

    public void collide(NewPlate plate, int globalX, int globalY, short crust) {
        int thisLocalX = this.getLocalX(globalX);
        int thisLocalY = this.getLocalY(globalY);
        int plateLocalX = plate.getLocalX(globalX);
        int plateLocalY = plate.getLocalY(globalY);
        double nx = thisLocalX - this.centerOfMassX - plateLocalX + plate.centerOfMassX;
        double ny = thisLocalY - this.centerOfMassY - plateLocalY + plate.centerOfMassY;
        double lenSqr = nx * nx + ny * ny;
        if (lenSqr <= 0) {
            return;
        }
        double length = Math.sqrt(lenSqr);
        nx /= (float) length;
        ny /= (float) length;
        final double dot = (this.velX - plate.velX) * nx + (this.velY - plate.velY) * ny;
        if (dot <= 0) {
            return;
        }
        double denom = lenSqr * (1.0 / this.mass + 1.0 / crust);
        double impulse = -dot / denom;
        this.accX += nx * impulse / this.mass;
        this.accY += ny * impulse / this.mass;
        plate.accX -= nx * impulse / (crust + plate.mass);
        plate.accY -= ny * impulse / (crust + plate.mass);
    }

    public boolean contains(int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        if (localX == Integer.MIN_VALUE) {
            return false;
        }
        int localY = this.getLocalY(globalY);
        if (localY == Integer.MIN_VALUE) {
            return false;
        }
        int index = this.getIndex(localX, localY);
        return this.heightmap[index] != Short.MIN_VALUE;
    }

    private int createSegment(int localX, int localY) {
        int locationIndex = localY * this.width + localX;
        int id = this.segmentData.size();
        if (this.segmentIdMap[locationIndex] < id) {
            return this.segmentIdMap[locationIndex];
        }
        boolean canGoLeft = localX > 0 && this.heightmap[locationIndex - 1] >= MainTectonics.CONTINENTAL_SHELF;
        boolean canGoRight = localX < this.width - 1 && this.heightmap[locationIndex + 1] >= MainTectonics.CONTINENTAL_SHELF;
        boolean canGoUp = localY > 0 && this.heightmap[locationIndex - this.width] >= MainTectonics.CONTINENTAL_SHELF;
        boolean canGoDown = localY < this.height - 1 && this.heightmap[locationIndex + this.width] >= MainTectonics.CONTINENTAL_SHELF;
        int neighbourId = id;
        // This point belongs to no segment yet.
        // However it might be a neighbour to some segment created earlier.
        // If such neighbour is found, associate this point with it.
        if (canGoLeft && this.segmentIdMap[locationIndex - 1] < id) {
            neighbourId = this.segmentIdMap[locationIndex - 1];
        }
        else if (canGoRight && this.segmentIdMap[locationIndex + 1] < id) {
            neighbourId = this.segmentIdMap[locationIndex + 1];
        }
        else if (canGoUp && this.segmentIdMap[locationIndex - this.width] < id) {
            neighbourId = this.segmentIdMap[locationIndex - this.width];
        }
        else if (canGoDown && this.segmentIdMap[locationIndex + this.width] < id) {
            neighbourId = this.segmentIdMap[locationIndex + this.width];
        }
        if (neighbourId < id) {
            this.segmentIdMap[locationIndex] = neighbourId;
            SegmentData segmentData = this.segmentData.get(neighbourId);
            ++segmentData.area;
            if (localY < segmentData.y0) {
                segmentData.y0 = localY;
            }
            if (localY > segmentData.y1) {
                segmentData.y1 = localY;
            }
            if (localX < segmentData.x0) {
                segmentData.x0 = localX;
            }
            if (localX > segmentData.x1) {
                segmentData.x1 = localX;
            }
            return neighbourId;
        }
        int linesProcessed;
        SegmentData data = new SegmentData(localX, localY, localX, localY, 0);
        IntList[] spans_todo = Noise.fill(new IntList[this.height], IntArrayList::new);
        IntList[] spans_done = Noise.fill(new IntList[this.height], IntArrayList::new);
        this.segmentIdMap[locationIndex] = id;
        spans_todo[localY].add(localX);
        spans_todo[localY].add(localX);
        do {
            linesProcessed = 0;
            for (int line = 0; line < this.height; ++line) {
                int start;
                int end;
                if (spans_todo[line].isEmpty()) {
                    continue;
                }
                // Find an unscanned span on this line.
                do {
                    end = spans_todo[line].removeInt(spans_todo[line].size() - 1);
                    start = spans_todo[line].removeInt(spans_todo[line].size() - 1);
                    // Reduce any done spans from this span.
                    for (int j = 0; j < spans_done[line].size(); j += 2) {
                        // Saved coordinates are AT the point
                        // that was included last to the span.
                        // That's why equalities matter.
                        if (start >= spans_done[line].get(j) &&
                            start <= spans_done[line].get(j + 1)) {
                            start = spans_done[line].get(j + 1) + 1;
                        }
                        if (end >= spans_done[line].get(j) &&
                            end <= spans_done[line].get(j + 1)) {
                            end = spans_done[line].get(j) - 1;
                        }
                    }
                    // Unsigned-ness hacking!
                    // Required to fix the underflow of end - 1.
                    start |= -(end >= this.width ? 1 : 0);
                    end -= end >= this.width ? 1 : 0;

                } while (start > end && !spans_todo[line].isEmpty());
                // Nothing to do here anymore...
                if (start > end) {
                    continue;
                }
                // Calculate line indices. Allow wrapping around map edges.
                final int row_above = line - 1 & -(line > 0 ? 1 : 0) | this.height - 1 & -(line == 0 ? 1 : 0);
                final int row_below = line + 1 & -(line < this.height - 1 ? 1 : 0);
                final int line_here = line * this.width;
                final int line_above = row_above * this.width;
                final int line_below = row_below * this.width;
                // Extend the beginning of line.
                while (start > 0 && this.segmentIdMap[line_here + start - 1] > id && this.heightmap[line_here + start - 1] >= MainTectonics.CONTINENTAL_SHELF) {
                    --start;
                    this.segmentIdMap[line_here + start] = id;
                    // Count volume of pixel...
                }
                // Extend the end of line.
                while (end < this.width - 1 && this.segmentIdMap[line_here + end + 1] > id && this.heightmap[line_here + end + 1] >= MainTectonics.CONTINENTAL_SHELF) {
                    ++end;
                    this.segmentIdMap[line_here + end] = id;
                    // Count volume of pixel...
                }
                // Check if should wrap around left edge.
                if (this.width == MainTectonics.WORLD_SIZE && start == 0 && this.segmentIdMap[line_here + this.width - 1] > id && this.heightmap[line_here + this.width - 1] >= MainTectonics.CONTINENTAL_SHELF) {
                    this.segmentIdMap[line_here + this.width - 1] = id;
                    spans_todo[line].add(this.width - 1);
                    spans_todo[line].add(this.width - 1);
                    // Count volume of pixel...
                }
                // Check if should wrap around right edge.
                if (this.width == MainTectonics.WORLD_SIZE && end == this.width - 1 && this.segmentIdMap[line_here] > id && this.heightmap[line_here] >= MainTectonics.CONTINENTAL_SHELF) {
                    this.segmentIdMap[line_here] = id;
                    spans_todo[line].add(0);
                    spans_todo[line].add(0);
                    // Count volume of pixel...
                }
                data.area += 1 + end - start; // Update segment area counter.
                // Record any changes in extreme dimensions.
                if (line < data.y0) {
                    data.y0 = line;
                }
                if (line > data.y1) {
                    data.y1 = line;
                }
                if (start < data.x0) {
                    data.x0 = start;
                }
                if (end > data.x1) {
                    data.x1 = end;
                }
                if (line > 0 || this.height == MainTectonics.WORLD_SIZE) {
                    for (int j = start; j <= end; ++j) {
                        if (this.segmentIdMap[line_above + j] > id && this.heightmap[line_above + j] >= MainTectonics.CONTINENTAL_SHELF) {
                            int a = j;
                            this.segmentIdMap[line_above + a] = id;
                            // Count volume of pixel...
                            while (++j < this.width && this.segmentIdMap[line_above + j] > id && this.heightmap[line_above + j] >= MainTectonics.CONTINENTAL_SHELF) {
                                this.segmentIdMap[line_above + j] = id;
                                // Count volume of pixel...
                            }
                            int b = --j; // Last point is invalid.
                            spans_todo[row_above].add(a);
                            spans_todo[row_above].add(b);
                            ++j; // Skip the last scanned point.
                        }
                    }
                }
                if (line < this.height - 1 || this.height == MainTectonics.WORLD_SIZE) {
                    for (int j = start; j <= end; ++j) {
                        if (this.segmentIdMap[line_below + j] > id && this.heightmap[line_below + j] >= MainTectonics.CONTINENTAL_SHELF) {
                            int a = j;
                            this.segmentIdMap[line_below + a] = id;
                            // Count volume of pixel...
                            while (++j < this.width && this.segmentIdMap[line_below + j] > id && this.heightmap[line_below + j] >= MainTectonics.CONTINENTAL_SHELF) {
                                this.segmentIdMap[line_below + j] = id;
                                // Count volume of pixel...
                            }
                            int b = --j; // Last point is invalid.
                            spans_todo[row_below].add(a);
                            spans_todo[row_below].add(b);
                            ++j; // Skip the last scanned point.
                        }
                    }
                }
                spans_done[line].add(start);
                spans_done[line].add(end);
                ++linesProcessed;
            }
        } while (linesProcessed > 0);
        this.segmentData.add(data);
        return id;
    }

    public void extendPlate(int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        if (localX != Integer.MIN_VALUE && localY != Integer.MIN_VALUE) {
            System.out.println("No need to extend plate");
            return;
        }
        double newX;
        int newWidth;
        if (localX == Integer.MIN_VALUE) {
            int dx1 = Math.abs((int) this.x - globalX - 1);
            int dx2 = MainTectonics.WORLD_SIZE - dx1;
            int dx = Math.max(dx1, dx2);
            int x1 = (int) this.x + this.width & MainTectonics.WORLD_SIZE - 1;
            int dw = x1 - globalX;
            dw = dw > 0 ? dw : MainTectonics.WORLD_SIZE + dw;
            newWidth = Math.min(dx, dw);
            if (dx == newWidth) {
                newX = this.x;
            }
            else {
                newX = x1 - dw & MainTectonics.WORLD_SIZE - 1;
            }
        }
        else {
            newX = this.x;
            newWidth = this.width;
        }
        double newY;
        int newHeight;
        if (localY == Integer.MIN_VALUE) {
            int dy1 = Math.abs((int) this.y - globalY - 1);
            int dy2 = MainTectonics.WORLD_SIZE - dy1;
            int dy = Math.max(dy1, dy2);
            int y1 = (int) this.y + this.height & MainTectonics.WORLD_SIZE - 1;
            int dh = y1 - globalY;
            dh = dh > 0 ? dh : MainTectonics.WORLD_SIZE + dh;
            newHeight = Math.min(dy, dh);
            if (dy == newHeight) {
                newY = this.y;
            }
            else {
                newY = y1 - dh & MainTectonics.WORLD_SIZE - 1;
            }
        }
        else {
            newY = this.y;
            newHeight = this.height;
        }
        int dx = (int) (this.x - newX) & MainTectonics.WORLD_SIZE - 1;
        int dy = (int) (this.y - newY) & MainTectonics.WORLD_SIZE - 1;
        if (newWidth != this.width) {
            if (newHeight != this.height) {
                //Needs to update both width and height
                short[] newHeightmap = new short[newWidth * newHeight];
                int[] newSegmentIdMap = new int[newWidth * newHeight];
                Arrays.fill(newHeightmap, Short.MIN_VALUE);
                Arrays.fill(newSegmentIdMap, Integer.MAX_VALUE);
                int yOffset = dy * newWidth;
                for (int i = 0; i < this.height; ++i) {
                    int from = i * this.width;
                    int to = i * newWidth + dx + yOffset;
                    System.arraycopy(this.heightmap, from, newHeightmap, to, this.width);
                    System.arraycopy(this.segmentIdMap, from, newSegmentIdMap, to, this.width);
                }
                if (this.x != newX || this.y != newY) {
                    for (int i = 0, len = this.segmentData.size(); i < len; ++i) {
                        SegmentData segmentData = this.segmentData.get(i);
                        segmentData.x0 += dx;
                        segmentData.x1 += dx;
                        segmentData.x0 &= MainTectonics.WORLD_SIZE - 1;
                        segmentData.x1 &= MainTectonics.WORLD_SIZE - 1;
                        segmentData.y0 += dy;
                        segmentData.y1 += dy;
                        segmentData.y0 &= MainTectonics.WORLD_SIZE - 1;
                        segmentData.y1 &= MainTectonics.WORLD_SIZE - 1;
                    }
                }
                this.x = newX;
                this.y = newY;
                this.width = newWidth;
                this.height = newHeight;
                this.heightmap = newHeightmap;
                this.segmentIdMap = newSegmentIdMap;
                assert this.getLocalX(globalX) != Integer.MIN_VALUE;
                assert this.getLocalY(globalY) != Integer.MIN_VALUE;
            }
            else {
                //Needs to update only width
                short[] newHeightmap = new short[newWidth * newHeight];
                int[] newSegmentIdMap = new int[newWidth * newHeight];
                Arrays.fill(newHeightmap, Short.MIN_VALUE);
                Arrays.fill(newSegmentIdMap, Integer.MAX_VALUE);
                for (int i = 0; i < this.height; ++i) {
                    int from = i * this.width;
                    int to = i * newWidth + dx;
                    System.arraycopy(this.heightmap, from, newHeightmap, to, this.width);
                    System.arraycopy(this.segmentIdMap, from, newSegmentIdMap, to, this.width);
                }
                if (this.x != newX) {
                    for (int i = 0, len = this.segmentData.size(); i < len; ++i) {
                        SegmentData segmentData = this.segmentData.get(i);
                        segmentData.x0 += dx;
                        segmentData.x1 += dx;
                        segmentData.x0 &= MainTectonics.WORLD_SIZE - 1;
                        segmentData.x1 &= MainTectonics.WORLD_SIZE - 1;
                    }
                }
                this.x = newX;
                this.width = newWidth;
                this.heightmap = newHeightmap;
                this.segmentIdMap = newSegmentIdMap;
                assert this.getLocalX(globalX) != Integer.MIN_VALUE;
            }
            return;
        }
        //Needs to update only height
        short[] newHeightmap = new short[newWidth * newHeight];
        int[] newSegmentIdMap = new int[newWidth * newHeight];
        Arrays.fill(newHeightmap, Short.MIN_VALUE);
        Arrays.fill(newSegmentIdMap, Integer.MAX_VALUE);
        int offset = dy * this.width;
        System.arraycopy(this.heightmap, 0, newHeightmap, offset, this.width * this.height);
        System.arraycopy(this.segmentIdMap, 0, newSegmentIdMap, offset, this.width * this.height);
        if (this.y != newY) {
            for (int i = 0, len = this.segmentData.size(); i < len; ++i) {
                SegmentData segmentData = this.segmentData.get(i);
                segmentData.y0 += dy;
                segmentData.y1 += dy;
                segmentData.y0 &= MainTectonics.WORLD_SIZE - 1;
                segmentData.y1 &= MainTectonics.WORLD_SIZE - 1;
            }
        }
        this.y = newY;
        this.height = newHeight;
        this.heightmap = newHeightmap;
        this.segmentIdMap = newSegmentIdMap;
        assert this.getLocalY(globalY) != Integer.MIN_VALUE;
//        int left = (int) this.x;
//        int top = (int) this.y;
//        int right = left + this.width - 1;
//        int bottom = top + this.height - 1;
//        globalX &= MainTectonics.WORLD_SIZE - 1;
//        globalY &= MainTectonics.WORLD_SIZE - 1;
//        // Calculate distance of new point from plate edges.
//        int left0 = left - globalX;
//        int right0 = (MainTectonics.WORLD_SIZE & -(globalX < left ? 1 : 0)) + globalX - right;
//        int top0 = top - globalY;
//        int bottom0 = (MainTectonics.WORLD_SIZE & -(globalY < top ? 1 : 0)) + globalY - bottom;
//        if (left0 < 0) {
//            left0 += Short.MAX_VALUE;
//        }
//        if (right0 < 0) {
//            right0 += Short.MAX_VALUE;
//        }
//        if (top0 < 0) {
//            top0 += Short.MAX_VALUE;
//        }
//        if (bottom0 < 0) {
//            bottom0 += Short.MAX_VALUE;
//        }
//        // Set larger of horizontal/vertical distance to zero.
//        // A valid distance is NEVER larger than world's side's length!
//        int d_lft = left0 & -(left0 < right0 ? 1 : 0) & -(left0 < MainTectonics.WORLD_SIZE ? 1 : 0);
//        int d_rgt = right0 & -(right0 <= left0 ? 1 : 0) & -(right0 < MainTectonics.WORLD_SIZE ? 1 : 0);
//        int d_top = top0 & -(top0 < bottom0 ? 1 : 0) & -(top0 < MainTectonics.WORLD_SIZE ? 1 : 0);
//        int d_btm = bottom0 & -(bottom0 <= top0 ? 1 : 0) & -(bottom0 < MainTectonics.WORLD_SIZE ? 1 : 0);
//        // Scale all changes to multiple of 8.
//        d_lft = (d_lft > 0 ? 1 : 0) + (d_lft >> 3) << 3;
//        d_rgt = (d_rgt > 0 ? 1 : 0) + (d_rgt >> 3) << 3;
//        d_top = (d_top > 0 ? 1 : 0) + (d_top >> 3) << 3;
//        d_btm = (d_btm > 0 ? 1 : 0) + (d_btm >> 3) << 3;
//        // Make sure plate doesn't grow bigger than the system it's in!
//        if (this.width + d_lft + d_rgt > MainTectonics.WORLD_SIZE) {
//            d_lft = 0;
//            d_rgt = MainTectonics.WORLD_SIZE - this.width;
//        }
//        if (this.height + d_top + d_btm > MainTectonics.WORLD_SIZE) {
//            d_top = 0;
//            d_btm = MainTectonics.WORLD_SIZE - this.height;
//        }
//        int oldWidth = this.width;
//        int oldHeight = this.height;
//        this.x -= d_lft;
//        this.x += this.x >= 0 ? 0 : MainTectonics.WORLD_SIZE;
//        this.width += d_lft + d_rgt;
//        this.y -= d_top;
//        this.y += this.y >= 0 ? 0 : MainTectonics.WORLD_SIZE;
//        this.height += d_top + d_btm;
//        short[] newHeightmap = new short[this.width * this.height];
//        int[] newSegmentIdMap = new int[this.width * this.height];
//        Arrays.fill(newHeightmap, Short.MIN_VALUE);
//        Arrays.fill(newSegmentIdMap, Integer.MAX_VALUE);
//        // copy old plate into new.
//        for (int j = 0; j < oldHeight; ++j) {
//            final int dstI = (d_top + j) * this.width + d_lft;
//            final int srcI = j * oldWidth;
//            System.arraycopy(this.heightmap, srcI, newHeightmap, dstI, oldWidth);
//            System.arraycopy(this.segmentIdMap, srcI, newSegmentIdMap, dstI, oldWidth);
//        }
//        this.heightmap = newHeightmap;
//        this.segmentIdMap = newSegmentIdMap;
//        // Shift all segment data to match new coordinates.
//        for (int s = 0; s < this.segmentData.size(); ++s) {
//            SegmentData segmentData = this.segmentData.get(s);
//            segmentData.x0 += d_lft;
//            segmentData.x1 += d_lft;
//            segmentData.y0 += d_top;
//            segmentData.y1 += d_top;
//        }
    }

    public SegmentData getCollisionInfo(int globalX, int globalY) {
        int index = this.getIndex(this.getLocalX(globalX), this.getLocalY(globalY));
        int segmentId = this.segmentIdMap[index];
        return this.segmentData.get(segmentId);
    }

    private short getCrustAmount(int globalX, int globalY) {
        int index = this.getIndex(globalX, globalY);
        if (index < 0 || index >= this.width * this.height) {
            return Short.MIN_VALUE;
        }
        return this.heightmap[index];
    }

    public int getHeight() {
        return this.height;
    }

    public short[] getHeightmap() {
        return this.heightmap;
    }

    private int getIndex(int localX, int localY) {
        if (MainTectonics.DEBUG) {
            if (localX == Integer.MIN_VALUE || localY == Integer.MIN_VALUE) {
                throw new IllegalStateException("Index is not valid!");
            }
        }
        return localY * this.width + localX;
    }

    /**
     * @param globalX The X coordinate on the global map. Might be unbounded.
     * @return The local X coordinate if within range, or {@link Integer#MIN_VALUE} if out of range.
     */
    private int getLocalX(int globalX) {
        int x0 = (int) this.x;
        int localX = globalX & MainTectonics.WORLD_SIZE - 1;
        if (localX < x0) {
            localX += MainTectonics.WORLD_SIZE;
        }
        localX -= x0;
        if (localX < 0) {
            throw new IllegalStateException("Zero");
        }
        if (localX >= this.width) {
            return Integer.MIN_VALUE;
        }
        return localX;
    }

    private int getLocalY(int globalY) {
        int y0 = (int) this.y;
        int localY = globalY & MainTectonics.WORLD_SIZE - 1;
        if (localY < y0) {
            localY += MainTectonics.WORLD_SIZE;
        }
        localY -= y0;
        if (localY < 0) {
            throw new IllegalStateException("Zero");
        }
        if (localY >= this.height) {
            return Integer.MIN_VALUE;
        }
        return localY;
    }

    public double getMass() {
        return this.mass;
    }

    public double getVelX() {
        return this.velX;
    }

    public double getVelY() {
        return this.velY;
    }

    public int getWidth() {
        return this.width;
    }

    public double getX() {
        return this.x;
    }

    public double getY() {
        return this.y;
    }

    public void makeSubdutionRamp(int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        short baseHeight = this.heightmap[this.getIndex(localX, localY)];
        //Up
        int index = this.getIndex(localX, localY - 1);
        if (0 <= index && index < this.width * this.height) {
            short localHeight = this.heightmap[index];
            if (localHeight > baseHeight && localHeight < MainTectonics.CONTINENTAL_SHELF) {
                this.heightmap[index] = (short) ((baseHeight + localHeight) / 2);
            }
        }
        //Down
        index = this.getIndex(localX, localY + 1);
        if (0 <= index && index < this.width * this.height) {
            short localHeight = this.heightmap[index];
            if (localHeight > baseHeight && localHeight < MainTectonics.CONTINENTAL_SHELF) {
                this.heightmap[index] = (short) ((baseHeight + localHeight) / 2);
            }
        }
        //Left
        index = this.getIndex(localX - 1, localY);
        if (0 <= index && index < this.width * this.height) {
            short localHeight = this.heightmap[index];
            if (localHeight > baseHeight && localHeight < MainTectonics.CONTINENTAL_SHELF) {
                this.heightmap[index] = (short) ((baseHeight + localHeight) / 2);
            }
        }
        //Right
        index = this.getIndex(localX + 1, localY);
        if (0 <= index && index < this.width * this.height) {
            short localHeight = this.heightmap[index];
            if (localHeight > baseHeight && localHeight < MainTectonics.CONTINENTAL_SHELF) {
                this.heightmap[index] = (short) ((baseHeight + localHeight) / 2);
            }
        }
    }

    public void move() {
        this.velX += this.accX;
        this.velY += this.accY;
        this.accX = 0;
        this.accY = 0;
        this.speed = Math.sqrt(this.velX * this.velX + this.velY * this.velY);
        this.velX /= this.speed;
        this.velY /= this.speed;
        double alpha = this.rotationDir * this.speed / (MainTectonics.WORLD_SIZE * 0.33);
        double cos = Math.cos(alpha * this.speed);
        double sin = Math.sin(alpha * this.speed);
        double newVx = this.velX * cos - this.velY * sin;
        double newVy = this.velY * cos + this.velX * sin;
        this.velX = newVx;
        this.velY = newVy;
        this.x += this.velX * this.speed;
        this.y += this.velY * this.speed;
        if (this.x < 0) {
            this.x += MainTectonics.WORLD_SIZE;
        }
        else if (this.x >= MainTectonics.WORLD_SIZE) {
            this.x -= MainTectonics.WORLD_SIZE;
        }
        if (this.y < 0) {
            this.y += MainTectonics.WORLD_SIZE;
        }
        else if (this.y >= MainTectonics.WORLD_SIZE) {
            this.y -= MainTectonics.WORLD_SIZE;
        }
    }

    public void resetSegments() {
        Arrays.fill(this.segmentIdMap, Integer.MAX_VALUE);
        this.segmentData.clear();
    }

    public void save(int id) throws IOException {
        BufferedImage image = new BufferedImage(MainTectonics.WORLD_SIZE, MainTectonics.WORLD_SIZE, BufferedImage.TYPE_INT_ARGB);
        for (int dy = 0, index = 0; dy < this.height; ++dy) {
            for (int dx = 0; dx < this.width; ++dx, ++index) {
                image.setRGB((int) (dx + this.x) & MainTectonics.WORLD_SIZE - 1, (int) (dy + this.y) & MainTectonics.WORLD_SIZE - 1, MainTectonics.getHeightmapColor(this.heightmap[index]));
            }
        }
        this.attachDebugInfo(image);
        ImageIO.write(image, "png", new File(String.format(Location.ROOT_FOLDER + "\\plate#%03d_%03d.png", id, this.imageId++)));
//        System.out.printf("Plate %3d has mass = %s\n", id, Metric.format(this.mass, 3, "mu"));
    }

    private void selectCollisionSegment(int globalX, int globalY) {
        int index = this.getIndex(this.getLocalX(globalX), this.getLocalY(globalY));
        this.activeContinent = this.segmentIdMap[index];
    }

    public void setCrust(int globalX, int globalY, short crust) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        if (localX == Integer.MIN_VALUE || localY == Integer.MIN_VALUE) {
            this.extendPlate(globalX, globalY);
        }
        int index = this.getIndex(this.getLocalX(globalX), this.getLocalY(globalY));
        this.mass -= this.heightmap[index];
        if (crust < MainTectonics.TRENCH_DEPTH) {
            this.heightmap[index] = Short.MIN_VALUE;
        }
        else {
            this.heightmap[index] = crust;
            this.mass += crust;
        }
    }

    public static class SegmentData {

        private int area;
        private int collisionCount;
        private int x0;
        private int x1;
        private int y0;
        private int y1;

        public SegmentData(int x0, int y0, int x1, int y1, int area) {
            this.x0 = x0;
            this.x1 = x1;
            this.y0 = y0;
            this.y1 = y1;
            this.area = area;
            this.collisionCount = 0;
        }

        public int getCollisionCount() {
            return this.collisionCount;
        }

        public float getCollisionRatio() {
            if (this.collisionCount == 0) {
                return 0;
            }
            return (float) this.collisionCount / this.area;
        }
    }
}
