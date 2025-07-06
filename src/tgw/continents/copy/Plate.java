package tgw.continents.copy;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.util.Arrays;

public class Plate {

    /**
     * Height limit that separates seas from dry land.
     */
    private static final float CONT_BASE = 1.0f;
    private static final int INITIAL_SPEED_X = 1;
    private static final int DEFORMATION_WEIGHT = 2;
    private static final boolean DEBUG = true;
    /**
     * Segment ID of the continent that's processed
     */
    private int activeContinent;
    /**
     * Bitmap of plate's soil's age, where age is the timestamp of creation
     */
    private int[] age;
    /**
     * X component of plate's centre of mass
     */
    private float centerOfMassX;
    /**
     * Y component of plate's centre of mass
     */
    private float centerOfMassY;
    /**
     * X component of plate's acceleration vector
     */
    private double dx;
    /**
     * Y component of plate's acceleration vector
     */
    private double dy;
    /**
     * Heightmap's height
     */
    private int height;
    /**
     * Bitmap of plate's structure/height
     */
    private float[] heightmap;
    private float mass;
    /**
     * Direction of plate's rotation, where 1 = CCW, -1 = CW
     */
    private float rot_dir;
    private final ObjectList<SegmentData> segData;
    /**
     * Segment ID of each piece of continental crust
     */
    private int[] segment;
    /**
     * X component of plate's direction unit vector
     */
    private double velX;
    /**
     * Y component of plate's direction unit vector
     */
    private double velY;
    private double velocity;
    /**
     * Heightmap's width
     */
    private int width;
    private final int worldSize;
    private double x;
    private double y;

    /**
     * Initializes plate with the supplied height map.
     *
     * @param m          Heightmap of terrain
     * @param w          Width of heightmap in pixels
     * @param h          Height of heightmap in pixels
     * @param _x         x coordinate of heightmap's top-left corner on world map
     * @param _y         y coordinate of heightmap's top-left corner on world map
     * @param world_size Length of world map in pixels (map is square)
     */
    public Plate(float[] m, int w, int h, int _x, int _y, int plate_age, int world_size) {
        this.segData = new ObjectArrayList<>();
        this.width = w;
        this.height = h;
        this.worldSize = world_size;
        this.mass = 0;
        this.x = _x;
        this.y = _y;
        this.centerOfMassX = 0;
        this.centerOfMassY = 0;
        this.dx = 0;
        this.dy = 0;
        final int A = w * h;
        final double angle = 2 * Math.PI * DumbMain.RANDOM.nextDouble();
        int k;
        if (m == null) {
            return;
        }
        this.heightmap = new float[A];
        this.age = new int[A];
        this.segment = new int[A];
        this.velocity = 1;
        this.rot_dir = DumbMain.RANDOM.nextBoolean() ? 1 : -1;
        this.velX = Math.cos(angle) * INITIAL_SPEED_X;
        this.velY = Math.sin(angle) * INITIAL_SPEED_X;
        Arrays.fill(this.segment, Integer.MAX_VALUE);
        for (int j = k = 0; j < this.height; ++j) {
            for (int i = 0; i < this.width; ++i, ++k) {
                // Clone map data and count crust mass.
                this.mass += this.heightmap[k] = m[k];
                // Calculate center coordinates weighted by mass.
                this.centerOfMassX += i * m[k];
                this.centerOfMassY += j * m[k];
                // Set the age of ALL points in this plate to same
                // value. The right thing to do would be to simulate
                // the generation of new oceanic crust as if the plate
                // had been moving to its current direction until all
                // plate's (oceanic) crust receive an age.
                this.age[k] = plate_age & -(m[k] > 0 ? 1 : 0);
            }
        }
        // Normalize center of mass coordinates.
        this.centerOfMassX /= this.mass;
        this.centerOfMassY /= this.mass;
    }

    /**
     * Increment collision counter of the continent at given location.
     *
     * @param globalX x coordinate of collision point on world map
     * @param globalY y coordinate of collision point on world map
     * @return Surface area of the collided continent (hack!)
     */
    public int addCollision(int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        if (DEBUG) {
            if (index >= this.width * this.height) {
                throw new IllegalStateException("Continental collision out of map bounds!");
            }
        }
        int seg = this.segment[index];
        if (seg >= this.segData.size()) {
            seg = this.createSegment(localX, localY);
        }
        if (DEBUG) {
            if (seg >= this.segData.size()) {
                throw new IllegalStateException("Could not create segment!");
            }
        }
        ++this.segData.get(seg).collisionCount;
        return this.segData.get(seg).area;
    }

    /**
     * Add crust to plate as result of continental collision.
     *
     * @param globalX X coordinate of new crust on global world map
     * @param globalY Y coordinate of new crust on global world map
     * @param crust   Amount of crust to add
     * @param time    Time of creation of new crust
     */
    public void addCrustByCollision(int globalX, int globalY, float crust, int time) {
        // Add crust. Extend plate if necessary.
        this.setCrust(globalX, globalY, this.getCrustAmount(globalX, globalY) + crust, time);
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        if (DEBUG) {
            if (index >= this.width * this.height) {
                System.out.println("Aggregation went overboard!");
                System.exit(1);
            }
        }
        this.segment[index] = this.activeContinent;
        SegmentData data = this.segData.get(this.activeContinent);
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

    /**
     * Simulates subduction of oceanic plate under this plate.
     * <p>
     * Subduction is simulated by calculating the distance on surface
     * that subducting sediment will travel under the plate until the
     * subducting slab has reached certain depth where the heat triggers
     * the melting and uprising of molten magma.
     *
     * @param globalX X coordinate of subduction on world map
     * @param globalY Y coordinate of subduction on world map
     * @param crust   Amount of sediment that subducts
     * @param time    Time of creation of new crust
     * @param velX    X direction of the subducting plate
     * @param velY    Y direction of the subducting plate
     */
    public void addCrustBySubduction(int globalX, int globalY, float crust, int time, double velX, double velY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        if (DEBUG) {
            // Should never be true!
            if (index >= this.width * this.height) {
                System.out.println("Subduction origin not on plate!");
                System.out.printf("%d, %d @ [%f, %f]x[%d, %d]\n", localX, localY, this.x, this.y, this.width, this.height);
                System.exit(1);
            }
        }
        // Take vector difference only between plates that move more or less
        // to same direction. This makes subduction direction behave better.
        if (this.velX * velX + this.velY * velY > 0) {
            velX -= this.velX;
            velY -= this.velY;
        }
        float offset = DumbMain.RANDOM.nextFloat();
        offset *= offset * offset * (2 * DumbMain.RANDOM.nextInt(2) - 1);
        velX = 10 * velX + 3 * offset;
        velY = 10 * velY + 3 * offset;
        localX += (int) velX;
        localY += (int) velY;
        if (this.width == this.worldSize) {
            localX &= this.width - 1;
        }
        if (this.height == this.worldSize) {
            localY &= this.height - 1;
        }
        index = localY * this.width + localX;
        if (index >= 0 && index < this.width * this.height && this.heightmap[index] > 0) {
            time = (int) ((this.heightmap[index] * this.age[index] + crust * time) / (this.heightmap[index] + crust));
            this.age[index] = time * (crust > 0 ? 1 : 0);
            this.heightmap[index] += crust;
            this.mass += crust;
        }
    }

    /**
     * Add continental crust from this plate as part of other plate.
     * <p>
     * Aggregation of two continents is the event where the collided
     * pieces of crust fuse together at the point of collision. It is
     * crucial to merge not only the collided pieces of crust but also
     * the entire continent that's part of the colliding tad of crust
     * However, because one plate can contain many islands and pieces of
     * continents, the merging must be done WITHOUT merging the entire
     * plate and all those continental pieces that have NOTHING to do with
     * the collision in question.
     *
     * @param plate   Reference to the receiving plate
     * @param globalX X coordinate of collision on world map
     * @param globalY Y coordinate of collision on world map
     * @return Amount of crust aggregated to destination plate
     */
    public float aggregateCrust(Plate plate, int globalX, int globalY) {
        int localX = this.getLocalX(globalX);
        int localY = this.getLocalY(globalY);
        int index = localY * this.width + localX;
        if (DEBUG) {
            if (index >= this.width * this.height) {
                System.out.println("Trying to aggregate beyond plate limits!");
                System.exit(1);
            }
        }
        int segmentId = this.getSegmentId(globalX, globalY, index);
        // One continent may have many points of collision. If one of them
        // causes continent to aggregate then all successive collisions and
        // attempts of aggregation would necessarily change nothing at all,
        // because the continent was removed from this plate earlier!
        SegmentData segmentData = this.segData.get(segmentId);
        if (segmentData.area == 0) {
            return 0;    // Do not process empty continents.
        }
        plate.selectCollisionSegment(globalX, globalY);
        // Wrap coordinates around world edges to safeguard subtractions.
        globalX += this.worldSize;
        globalY += this.worldSize;
        float oldMass = this.mass;
        // Add all the collided continent's crust to destination plate.
        for (int y = segmentData.y0; y <= segmentData.y1; ++y) {
            for (int x = segmentData.x0; x <= segmentData.x1; ++x) {
                int i = y * this.width + x;
                if (this.segment[i] == segmentId && this.heightmap[i] > 0) {
                    plate.addCrustByCollision(globalX + x - localX, globalY + y - localY, this.heightmap[i], this.age[i]);
                    this.mass -= this.heightmap[i];
                    this.heightmap[i] = 0;
                }
            }
        }
        segmentData.area = 0; // Mark segment as non-existent.
        return oldMass - this.mass;
    }

    /**
     * Decrease the speed of plate amount relative to its total mass.
     * <p>
     * Method decreases the speed of plate due to friction that occurs when
     * two plates collide. The amount of reduction depends on the amount
     * of mass that causes friction (i.e. that has collided) compared to
     * the total mass of the plate. Thus, big chunk of crust colliding to
     * a small plate will halt it but have little effect on a huge plate.
     *
     * @param deformed_mass Amount of mass deformed in the collision
     */
    public void applyFriction(float deformed_mass) {
        // Remove the energy that deformation consumed from plate's kinetic
        // energy: F - dF = ma - dF => a = dF/m.
        if (this.mass > 0) {
            double vel_dec = DEFORMATION_WEIGHT * deformed_mass / this.mass;
            vel_dec = Math.min(vel_dec, this.velocity);
            // Altering the source variable causes the order of calls to
            // this function to have difference when it shouldn't!
            // However, it's a hack well worth the outcome. :)
            this.velocity -= vel_dec;
        }
    }

    /**
     * Method collides two plates according to Newton's laws of motion.
     * <p>
     * The velocity and direction of both plates are updated using
     * impulse forces following the collision according to Newton's laws
     * of motion. Deformations are not applied but energy consumed by the
     * deformation process is taken away from plate's momentum.
     *
     * @param plate         NewPlate to test against
     * @param globalX       X coordinate of collision on world map
     * @param globalY       Y coordinate of collision on world map
     * @param collidingMass Amount of colliding mass from source plate
     */
    public void collide(Plate plate, int globalX, int globalY, float collidingMass) {
        // Calculate the normal to the curve/line at collision point.
        // The normal will point into plate B i.e. the "other" plate.
        //
        // Plates that wrap over world edges can mess the normal vector.
        // This could be solved by choosing the normal vector that points the
        // shortest path beween mass centers but this causes problems when
        // plates are like heavy metal balls at a long rod and one plate's ball
        // collides at the further end of other plate's rod. Sure, this is
        // nearly never occurring situation but if we can easily do better then
        // why not do it?
        //
        // Better way is to select that normal vector that points along the
        // line that passes nearest the point of collision. Because point's
        // distance from line segment is relatively cumbersome to perform, the
        // vector is constructed as the sum of vectors <massCenterA, P> and
        // <P, massCenterB>. This solution works because collisions always
        // happen in the overlapping region of the two plates.
        int thisLocalX = this.getLocalX(globalX);
        int thisLocalY = this.getLocalY(globalY);
        int plateLocalX = plate.getLocalX(globalX);
        int plateLocalY = plate.getLocalY(globalY);
        if (DEBUG) {
            int thisIndex = thisLocalY * this.width + thisLocalX;
            int plateIndex = plateLocalY * plate.width + plateLocalX;
            if (thisIndex >= this.width * this.height || plateIndex >= plate.width * plate.height) {
                System.out.printf("@%d, %d: out of colliding map's bounds!\n", globalX, globalY);
                System.exit(1);
            }
        }
        float nx = thisLocalX - this.centerOfMassX - plateLocalX + plate.centerOfMassX;
        float ny = thisLocalY - this.centerOfMassY - plateLocalY + plate.centerOfMassY;
        float lenSqr = nx * nx + ny * ny;
        if (lenSqr <= 0) {
            return; // Avoid division by zero!
        }
        // Scaling is required at last when impulses are added to plates!
        double length = Math.sqrt(lenSqr);
        nx /= (float) length;
        ny /= (float) length;
        // Compute relative velocity between plates at the collision point.
        // Because torque is not included, calc simplifies to v_ab = v_a - v_b.
        // Get the dot product of relative velocity vector and collision vector.
        // Then get the projection of v_ab along collision vector.
        // Note that vector n must be a unit vector!
        final double dot = (this.velX - plate.velX) * nx + (this.velY - plate.velY) * ny;
        if (dot <= 0) {
            return; // Exit if objects are moving away from each other.
        }
        // Calculate the denominator of impulse: n . n * (1 / m_1 + 1 / m_2).
        // Use the mass of the colliding crust for the "donator" plate.
        double denom = lenSqr * (1.0 / this.mass + 1.0 / collidingMass);
        // Calculate force of impulse.
        double impulse = -dot / denom;
        // Compute final change of trajectory.
        // The plate that is the "giver" of the impulse should receive a
        // force according to its pre-collision mass, not the current mass!
        this.dx += nx * impulse / this.mass;
        this.dy += ny * impulse / this.mass;
        plate.dx -= nx * impulse / (collidingMass + plate.mass);
        plate.dy -= ny * impulse / (collidingMass + plate.mass);
    }

    /**
     * Separate a continent at (X, Y) to its own partition.
     * <p>
     * Method analyzes the pixels 4-ways adjacent at the given location
     * and labels all connected continental points with same segment ID.
     *
     * @param localX X offset on the local heightmap
     * @param localY Y offset on the local heightmap
     * @return ID of created segment on success, otherwise -1
     */
    private int createSegment(int localX, int localY) {
        int locationIndex = localY * this.width + localX;
        int id = this.segData.size();
        if (this.segment[locationIndex] < id) {
            return this.segment[locationIndex];
        }
        boolean canGoLeft = localX > 0 && this.heightmap[locationIndex - 1] >= CONT_BASE;
        boolean canGoRight = localX < this.width - 1 && this.heightmap[locationIndex + 1] >= CONT_BASE;
        boolean canGoUp = localY > 0 && this.heightmap[locationIndex - this.width] >= CONT_BASE;
        boolean canGoDown = localY < this.height - 1 && this.heightmap[locationIndex + this.width] >= CONT_BASE;
        int neighbourId = id;
        // This point belongs to no segment yet.
        // However it might be a neighbour to some segment created earlier.
        // If such neighbour is found, associate this point with it.
        if (canGoLeft && this.segment[locationIndex - 1] < id) {
            neighbourId = this.segment[locationIndex - 1];
        }
        else if (canGoRight && this.segment[locationIndex + 1] < id) {
            neighbourId = this.segment[locationIndex + 1];
        }
        else if (canGoUp && this.segment[locationIndex - this.width] < id) {
            neighbourId = this.segment[locationIndex - this.width];
        }
        else if (canGoDown && this.segment[locationIndex + this.width] < id) {
            neighbourId = this.segment[locationIndex + this.width];
        }
        if (neighbourId < id) {
            this.segment[locationIndex] = neighbourId;
            ++this.segData.get(neighbourId).area;
            if (localY < this.segData.get(neighbourId).y0) {
                this.segData.get(neighbourId).y0 = localY;
            }
            if (localY > this.segData.get(neighbourId).y1) {
                this.segData.get(neighbourId).y1 = localY;
            }
            if (localX < this.segData.get(neighbourId).x0) {
                this.segData.get(neighbourId).x0 = localX;
            }
            if (localX > this.segData.get(neighbourId).x1) {
                this.segData.get(neighbourId).x1 = localX;
            }
            return neighbourId;
        }
        int lines_processed;
        SegmentData data = new SegmentData(localX, localY, localX, localY, 0);
        IntList[] spans_todo = Noise.fill(new IntList[this.height], IntArrayList::new);
        IntList[] spans_done = Noise.fill(new IntList[this.height], IntArrayList::new);
        this.segment[locationIndex] = id;
        spans_todo[localY].add(localX);
        spans_todo[localY].add(localX);
        do {
            lines_processed = 0;
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
                while (start > 0 && this.segment[line_here + start - 1] > id && this.heightmap[line_here + start - 1] >= CONT_BASE) {
                    --start;
                    this.segment[line_here + start] = id;
                    // Count volume of pixel...
                }
                // Extend the end of line.
                while (end < this.width - 1 && this.segment[line_here + end + 1] > id && this.heightmap[line_here + end + 1] >= CONT_BASE) {
                    ++end;
                    this.segment[line_here + end] = id;
                    // Count volume of pixel...
                }
                // Check if should wrap around left edge.
                if (this.width == this.worldSize && start == 0 &&
                    this.segment[line_here + this.width - 1] > id &&
                    this.heightmap[line_here + this.width - 1] >= CONT_BASE) {
                    this.segment[line_here + this.width - 1] = id;
                    spans_todo[line].add(this.width - 1);
                    spans_todo[line].add(this.width - 1);
                    // Count volume of pixel...
                }
                // Check if should wrap around right edge.
                if (this.width == this.worldSize && end == this.width - 1 &&
                    this.segment[line_here] > id &&
                    this.heightmap[line_here] >= CONT_BASE) {
                    this.segment[line_here] = id;
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
                if (line > 0 || this.height == this.worldSize) {
                    for (int j = start; j <= end; ++j) {
                        if (this.segment[line_above + j] > id &&
                            this.heightmap[line_above + j] >= CONT_BASE) {
                            int a = j;
                            this.segment[line_above + a] = id;
                            // Count volume of pixel...
                            while (++j < this.width && this.segment[line_above + j] > id && this.heightmap[line_above + j] >= CONT_BASE) {
                                this.segment[line_above + j] = id;
                                // Count volume of pixel...
                            }
                            int b = --j; // Last point is invalid.
                            spans_todo[row_above].add(a);
                            spans_todo[row_above].add(b);
                            ++j; // Skip the last scanned point.
                        }
                    }
                }
                if (line < this.height - 1 || this.height == this.worldSize) {
                    for (int j = start; j <= end; ++j) {
                        if (this.segment[line_below + j] > id &&
                            this.heightmap[line_below + j] >= CONT_BASE) {
                            int a = j;
                            this.segment[line_below + a] = id;
                            // Count volume of pixel...
                            while (++j < this.width && this.segment[line_below + j] > id && this.heightmap[line_below + j] >= CONT_BASE) {
                                this.segment[line_below + j] = id;
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
                ++lines_processed;
            }
        } while (lines_processed > 0);
        this.segData.add(data);
        return id;
    }

    /**
     * Apply plate wide erosion algorithm.
     * <p>
     * Plates total mass and the center of mass are updated.
     *
     * @param lower_bound Sets the limit below which there's no erosion
     */
    public void erode(float lower_bound) {
        IntList sources_data = new IntArrayList();
        IntList sinks_data = new IntArrayList();
        IntList sources = sources_data;
        IntList sinks = sinks_data;
        float[] tmp = new float[this.width * this.height];
        System.arraycopy(this.heightmap, 0, tmp, 0, this.width * this.height);
        // Find all tops.
        for (int y = 0; y < this.height; ++y) {
            for (int x = 0; x < this.width; ++x) {
                final int index = y * this.width + x;
                if (this.heightmap[index] < lower_bound) {
                    continue;
                }
                // Build masks for accessible directions (4-way).
                // Allow wrapping around map edges if plate has world wide dimensions.
                int w_mask = -((x > 0 ? 1 : 0) | (this.width == this.worldSize ? 1 : 0));
                int e_mask = -((x < this.width - 1 ? 1 : 0) | (this.width == this.worldSize ? 1 : 0));
                int n_mask = -((y > 0 ? 1 : 0) | (this.height == this.worldSize ? 1 : 0));
                int s_mask = -((y < this.height - 1 ? 1 : 0) | (this.height == this.worldSize ? 1 : 0));
                // Calculate the x and y offset of neighbour directions.
                // If neighbour is out of plate edges, set it to zero. This protects
                // map memory reads from segment faulting.
                int w = this.worldSize + x - 1 & this.worldSize - 1 & w_mask;
                int e = this.worldSize + x + 1 & this.worldSize - 1 & e_mask;
                int n = this.worldSize + y - 1 & this.worldSize - 1 & n_mask;
                int s = this.worldSize + y + 1 & this.worldSize - 1 & s_mask;
                // Calculate offsets within map memory.
                w = y * this.width + w;
                e = y * this.width + e;
                n = n * this.width + x;
                s = s * this.width + x;
                // Extract neighbours heights. Apply validity filtering: 0 is invalid.
                float w_crust = this.heightmap[w] * (w_mask & (this.heightmap[w] < this.heightmap[index] ? 1 : 0));
                float e_crust = this.heightmap[e] * (e_mask & (this.heightmap[e] < this.heightmap[index] ? 1 : 0));
                float n_crust = this.heightmap[n] * (n_mask & (this.heightmap[n] < this.heightmap[index] ? 1 : 0));
                float s_crust = this.heightmap[s] * (s_mask & (this.heightmap[s] < this.heightmap[index] ? 1 : 0));
                // This location is either at the edge of the plate or it is not the
                // tallest of its neightbours. Don't start a river from here.
                if (w_crust * e_crust * n_crust * s_crust == 0) {
                    continue;
                }
                sources.add(index);
            }
        }
        int[] isDone = new int[this.width * this.height];
        // From each top, start flowing water along the steepest slope.
        while (!sources.isEmpty()) {
            while (!sources.isEmpty()) {
                final int index = sources.getInt(sources.size() - 1);
                final int y = index / this.width;
                final int x = index - y * this.width;
                sources.removeInt(sources.size() - 1);
                if (this.heightmap[index] < lower_bound) {
                    continue;
                }
                // Build masks for accessible directions (4-way).
                // Allow wrapping around map edges if plate has worldwide dimensions.
                int w_mask = -((x > 0 ? 1 : 0) | (this.width == this.worldSize ? 1 : 0));
                int e_mask = -((x < this.width - 1 ? 1 : 0) | (this.width == this.worldSize ? 1 : 0));
                int n_mask = -((y > 0 ? 1 : 0) | (this.height == this.worldSize ? 1 : 0));
                int s_mask = -((y < this.height - 1 ? 1 : 0) | (this.height == this.worldSize ? 1 : 0));
                // Calculate the x and y offset of neighbour directions.
                // If neighbour is out of plate edges, set it to zero. This protects
                // map memory reads from segment faulting.
                int w = this.worldSize + x - 1 & this.worldSize - 1 & w_mask;
                int e = this.worldSize + x + 1 & this.worldSize - 1 & e_mask;
                int n = this.worldSize + y - 1 & this.worldSize - 1 & n_mask;
                int s = this.worldSize + y + 1 & this.worldSize - 1 & s_mask;
                // Calculate offsets within map memory.
                w = y * this.width + w;
                e = y * this.width + e;
                n = n * this.width + x;
                s = s * this.width + x;
                // Extract neighbours heights. Apply validity filtering: 0 is invalid.
                float w_crust = this.heightmap[w] * (w_mask & (this.heightmap[w] < this.heightmap[index] ? 1 : 0));
                float e_crust = this.heightmap[e] * (e_mask & (this.heightmap[e] < this.heightmap[index] ? 1 : 0));
                float n_crust = this.heightmap[n] * (n_mask & (this.heightmap[n] < this.heightmap[index] ? 1 : 0));
                float s_crust = this.heightmap[s] * (s_mask & (this.heightmap[s] < this.heightmap[index] ? 1 : 0));
                // If this is the lowest part of its neighbourhood, stop.
                if (w_crust + e_crust + n_crust + s_crust == 0) {
                    continue;
                }
                w_crust += (w_crust == 0 ? 1 : 0) * this.heightmap[index];
                e_crust += (e_crust == 0 ? 1 : 0) * this.heightmap[index];
                n_crust += (n_crust == 0 ? 1 : 0) * this.heightmap[index];
                s_crust += (s_crust == 0 ? 1 : 0) * this.heightmap[index];
                // Find the lowest neighbour.
                float lowest_crust = w_crust;
                int dest = index - 1;
                if (e_crust < lowest_crust) {
                    lowest_crust = e_crust;
                    dest = index + 1;
                }
                if (n_crust < lowest_crust) {
                    lowest_crust = n_crust;
                    dest = index - this.width;
                }
                if (s_crust < lowest_crust) {
                    lowest_crust = s_crust;
                    dest = index + this.width;
                }
                // if it's not handled yet, add it as new sink.
                if (dest < this.width * this.height && isDone[dest] == 0) {
                    sinks.add(dest);
                    isDone[dest] = 1;
                }
                // Erode this location with the water flow.
                tmp[index] -= (tmp[index] - lower_bound) * 0.2f;
            }
            IntList v_tmp = sources;
            sources = sinks;
            sinks = v_tmp;
            sinks.clear();
        }
        // Add random noise (10 %) to heightmap.
        for (int i = 0; i < this.width * this.height; ++i) {
            float alpha = 0.2f * DumbMain.RANDOM.nextFloat();
            tmp[i] += 0.1f * tmp[i] - alpha * tmp[i];
        }
        System.arraycopy(tmp, 0, this.heightmap, 0, this.width * this.height);
        Arrays.fill(tmp, 0);
        this.mass = 0;
        this.centerOfMassX = this.centerOfMassY = 0;
        for (int y = 0; y < this.height; ++y) {
            for (int x = 0; x < this.width; ++x) {
                final int index = y * this.width + x;
                this.mass += this.heightmap[index];
                tmp[index] += this.heightmap[index];
                // Update the center coordinates weighted by mass.
                this.centerOfMassX += x * this.heightmap[index];
                this.centerOfMassY += y * this.heightmap[index];
                if (this.heightmap[index] < lower_bound) {
                    continue;
                }
                // Build masks for accessible directions (4-way).
                // Allow wrapping around map edges if plate has worldwide dimensions.
                int w_mask = -((x > 0 ? 1 : 0) | (this.width == this.worldSize ? 1 : 0));
                int e_mask = -((x < this.width - 1 ? 1 : 0) | (this.width == this.worldSize ? 1 : 0));
                int n_mask = -((y > 0 ? 1 : 0) | (this.height == this.worldSize ? 1 : 0));
                int s_mask = -((y < this.height - 1 ? 1 : 0) | (this.height == this.worldSize ? 1 : 0));
                // Calculate the x and y offset of neighbour directions.
                // If neighbour is out of plate edges, set it to zero. This protects
                // map memory reads from segment faulting.
                int w = this.worldSize + x - 1 & this.worldSize - 1 & w_mask;
                int e = this.worldSize + x + 1 & this.worldSize - 1 & e_mask;
                int n = this.worldSize + y - 1 & this.worldSize - 1 & n_mask;
                int s = this.worldSize + y + 1 & this.worldSize - 1 & s_mask;
                // Calculate offsets within map memory.
                w = y * this.width + w;
                e = y * this.width + e;
                n = n * this.width + x;
                s = s * this.width + x;
                // Extract neighbours heights. Apply validity filtering: 0 is invalid.
                float w_crust = this.heightmap[w] * (w_mask & (this.heightmap[w] < this.heightmap[index] ? 1 : 0));
                float e_crust = this.heightmap[e] * (e_mask & (this.heightmap[e] < this.heightmap[index] ? 1 : 0));
                float n_crust = this.heightmap[n] * (n_mask & (this.heightmap[n] < this.heightmap[index] ? 1 : 0));
                float s_crust = this.heightmap[s] * (s_mask & (this.heightmap[s] < this.heightmap[index] ? 1 : 0));
                // This location has no neighbours (ARTIFACT!) or it is the lowest
                // part of its area. In either case the work here is done.
                if (w_crust + e_crust + n_crust + s_crust == 0) {
                    continue;
                }
                // The steeper the slope, the more water flows along it.
                // The more downhill (sources), the more water flows to here.
                // 1+1+10 = 12, avg = 4, stdev = sqrt((3*3+3*3+6*6)/3) = 4.2, var = 18,
                //	1*1+1*1+10*10 = 102, 102/4.2=24
                // 1+4+7 = 12, avg = 4, stdev = sqrt((3*3+0*0+3*3)/3) = 2.4, var = 6,
                //	1*1+4*4+7*7 = 66, 66/2.4 = 27
                // 4+4+4 = 12, avg = 4, stdev = sqrt((0*0+0*0+0*0)/3) = 0, var = 0,
                //	4*4+4*4+4*4 = 48, 48/0 = inf -> 48
                // If there's a source slope of height X then it will always cause
                // water erosion of amount Y. Then again from one spot only so much
                // water can flow.
                // Thus, the calculated non-linear flow value for this location is
                // multiplied by the "water erosion" constant.
                // The result is max(result, 1.0). New height of this location could
                // be e.g. h_lowest + (1 - 1 / result) * (h_0 - h_lowest).
                // Calculate the difference in height between this point and its
                // nbours that are lower than this point.
                float w_diff = this.heightmap[index] - w_crust;
                float e_diff = this.heightmap[index] - e_crust;
                float n_diff = this.heightmap[index] - n_crust;
                float s_diff = this.heightmap[index] - s_crust;
                float min_diff = w_diff;
                min_diff -= (min_diff - e_diff) * (e_diff < min_diff ? 1 : 0);
                min_diff -= (min_diff - n_diff) * (n_diff < min_diff ? 1 : 0);
                min_diff -= (min_diff - s_diff) * (s_diff < min_diff ? 1 : 0);
                // Calculate the sum of difference between lower neighbours and
                // the TALLEST lower neighbour.
                float diff_sum = (w_diff - min_diff) * (w_crust > 0 ? 1 : 0) +
                                 (e_diff - min_diff) * (e_crust > 0 ? 1 : 0) +
                                 (n_diff - min_diff) * (n_crust > 0 ? 1 : 0) +
                                 (s_diff - min_diff) * (s_crust > 0 ? 1 : 0);
                if (DEBUG) {
                    if (diff_sum < 0) {
                        System.out.println("Erosion differense sum is negative!");
                        System.out.printf("%f > %f %f %f %f\n", min_diff, w_diff, e_diff, n_diff, s_diff);
                        System.exit(1);
                    }
                }
                if (diff_sum < min_diff) {
                    // There's NOT enough room in neighbours to contain all the
                    // crust from this peak so that it would be as tall as its
                    // tallest lower neighbour. Thus first step is make ALL
                    // lower neighbours and this point equally tall.
                    tmp[w] += (w_diff - min_diff) * (w_crust > 0 ? 1 : 0);
                    tmp[e] += (e_diff - min_diff) * (e_crust > 0 ? 1 : 0);
                    tmp[n] += (n_diff - min_diff) * (n_crust > 0 ? 1 : 0);
                    tmp[s] += (s_diff - min_diff) * (s_crust > 0 ? 1 : 0);
                    tmp[index] -= min_diff;
                    min_diff -= diff_sum;
                    // Spread the remaining crust equally among all lower nbours.
                    min_diff /= 1 + (w_crust > 0 ? 1 : 0) + (e_crust > 0 ? 1 : 0) + (n_crust > 0 ? 1 : 0) + (s_crust > 0 ? 1 : 0);
                    tmp[w] += min_diff * (w_crust > 0 ? 1 : 0);
                    tmp[e] += min_diff * (e_crust > 0 ? 1 : 0);
                    tmp[n] += min_diff * (n_crust > 0 ? 1 : 0);
                    tmp[s] += min_diff * (s_crust > 0 ? 1 : 0);
                    tmp[index] += min_diff;
                }
                else {
                    float unit = min_diff / diff_sum;
                    // Remove all crust from this location making it as tall as
                    // its tallest lower neighbour.
                    tmp[index] -= min_diff;
                    // Spread all removed crust among all other lower neighbours.
                    tmp[w] += unit * (w_diff - min_diff) * (w_crust > 0 ? 1 : 0);
                    tmp[e] += unit * (e_diff - min_diff) * (e_crust > 0 ? 1 : 0);
                    tmp[n] += unit * (n_diff - min_diff) * (n_crust > 0 ? 1 : 0);
                    tmp[s] += unit * (s_diff - min_diff) * (s_crust > 0 ? 1 : 0);
                }
            }
        }
        this.heightmap = tmp;
        if (this.mass > 0) {
            this.centerOfMassX /= this.mass;
            this.centerOfMassY /= this.mass;
        }
    }

    public int[] getAgeMap() {
        return this.age;
    }

    /**
     * Retrieve collision statistics of continent at given location.
     *
     * @param globalX X coordinate of collision point on world map
     * @param globalY Y coordinate of collision point on world map
     * @param count   Destination for the count of collisions
     * @param ratio   Destination for the % of area the collided
     */
    public SegmentData getCollisionInfo(int globalX, int globalY) {
        int index = this.getIndex(globalX, globalY);
        int seg = this.segData.size();
        if (DEBUG) {
            if (index >= this.width * this.height) {
                System.out.println("getCollisionInfo: out of map bounds!");
                System.exit(1);
            }
        }
        seg = this.segment[index];
        if (DEBUG) {
            if (seg >= this.segData.size()) {
                System.out.println("getCollisionInfo: no segment found!");
                System.exit(1);
            }
        }
        return this.segData.get(seg);
    }

    /**
     * Retrieve the surface area of continent lying at desired location.
     *
     * @param globalX X coordinate of collision on world map
     * @param globalY Y coordinate of collision on world map
     * @return Area of continent at desired location
     */
    public int getContinentArea(int globalX, int globalY) {
        int index = this.getIndex(globalX, globalY);
        if (DEBUG) {
            if (index >= this.width * this.height) {
                System.out.println("getContinentArea: out of map bounds!");
                System.exit(1);
            }
            if (this.segment[index] >= this.segData.size()) {
                System.out.println("getContinentArea: no segment found!");
                System.exit(1);
            }
        }
        return this.segData.get(this.segment[index]).area;
    }

    /**
     * Get the timestamp of plate's crustal material at some location.
     *
     * @param globalX X offset on the global world map
     * @param globalY Y offset on the global world map
     * @return Timestamp of creation of crust at the location or 0 if no crust
     */
    public int getCrustAge(int globalX, int globalY) {
        int index = this.getIndex(globalX, globalY);
        if (index < 0 || index >= this.width * this.height) {
            return 0;
        }
        return this.age[index];
    }

    /**
     * Get the amount of plate's crustal material at some location.
     *
     * @param globalX X offset on the global world map
     * @param globalY Y offset on the global world map
     * @return Amount of crust at requested location
     */
    public float getCrustAmount(int globalX, int globalY) {
        int index = this.getIndex(globalX, globalY);
        if (index < 0 || index >= this.width * this.height) {
            return 0;
        }
        return this.heightmap[index];
    }

    public int getHeight() {
        return this.height;
    }

    public float[] getHeightmap() {
        return this.heightmap;
    }

    private int getIndex(int globalX, int globalY) {
        return this.getLocalY(globalY) * this.width + this.getLocalX(globalX);
    }

    private int getLocalX(int globalX) {
        int x0 = (int) this.x;
        int localX = globalX & this.worldSize - 1;
        if (localX < x0) {
            localX += this.worldSize;
        }
        localX -= x0;
        return localX;
    }

    private int getLocalY(int globalY) {
        int y0 = (int) this.y;
        int localY = globalY & this.worldSize - 1;
        if (localY < y0) {
            localY += this.worldSize;
        }
        localY -= y0;
        return localY;
    }

    public float getMass() {
        return this.mass;
    }

    public double getMomentum() {
        return this.mass * this.velocity;
    }

    private int getSegmentId(int globalX, int globalY, int index) {
        final int segmentId = this.segment[index];
        // This check forces the caller to do things in proper order!
        //
        // Usually continents collide at several locations simultaneously.
        // Thus, if this segment that is being merged now is removed from
        // segmentation bookkeeping, then the next point of collision that is
        // processed during the same iteration step would cause the test
        // below to be true and system would experience a premature abort.
        //
        // Therefore, segmentation bookkeeping is left intact. It doesn't
        // cause significant problems because all crust is cleared and empty
        // points are not processed at all.
        if (DEBUG) {
            if (segmentId >= this.segData.size()) {
                throw new IllegalStateException("Trying to aggregate without deforming first at [" + globalX + ", " + globalY + "]");
            }
        }
        return segmentId;
    }

    public double getVelX() {
        return this.velX;
    }

    public double getVelY() {
        return this.velY;
    }

    public double getVelocity() {
        return this.velocity;
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

    public boolean isEmpty() {
        return this.mass <= 0;
    }

    /**
     * Moves the plate along its trajectory
     */
    public void move() {
        // Apply any new impulses to the plate's trajectory.
        this.velX += this.dx;
        this.velY += this.dy;
        this.dx = 0;
        this.dy = 0;
        // Force direction of plate to be unit vector.
        // Update velocity so that the distance of movement doesn't change.
        this.velocity = Math.sqrt(this.velX * this.velX + this.velY * this.velY);
        this.velX /= this.velocity;
        this.velY /= this.velocity;
        // Apply some circular motion to the plate.
        // Force the radius of the circle to remain fixed by adjusting
        // angular velocity (which depends on plate's velocity).
        double alpha = this.rot_dir * this.velocity / (this.worldSize * 0.33);
        double cos = Math.cos(alpha * this.velocity);
        double sin = Math.sin(alpha * this.velocity);
        double newVx = this.velX * cos - this.velY * sin;
        double newVy = this.velY * cos + this.velX * sin;
        this.velX = newVx;
        this.velY = newVy;
        // Location modulations into range [0, world_side[ are a have to!
        // If left undone SOMETHING WILL BREAK DOWN SOMEWHERE in the code!
        if (DEBUG) {
            if (this.x < 0 || this.x > this.worldSize || this.y < 0 || this.y > this.worldSize) {
                System.out.println("Location coordinates out of world map bounds (PRE)!");
                System.exit(1);
            }
        }
        this.x += this.velX * this.velocity;
        if (this.x < 0) {
            this.x += this.worldSize;
        }
        else if (this.x >= this.worldSize) {
            this.x -= this.worldSize;
        }
        this.y += this.velY * this.velocity;
        if (this.y < 0) {
            this.y += this.worldSize;
        }
        else if (this.y >= this.worldSize) {
            this.y -= this.worldSize;
        }
        if (DEBUG) {
            if (this.x < 0 || this.x > this.worldSize || this.y < 0 || this.y > this.worldSize) {
                System.out.println("Location coordinates out of world map bounds (POST)!");
                System.out.printf("%f, %f, %f; %f, %f\n", this.velX, this.velY, this.velocity, this.x, this.y);
                System.exit(1);
            }
        }
    }

    /**
     * Clear any earlier continental crust partitions.
     * <p>
     * NewPlate has an internal bookkeeping of distinct areas of continental
     * crust for more realistic collision response. However, as the number
     * of collisions that plate experiences grows, so does the bookkeeping
     * of a continent become more and more inaccurate. Finally, it results
     * in striking artefacts that cannot be overlooked.
     * <p>
     * To alleviate this problem without the need of per iteration
     * recalculations plate supplies caller a method to reset its
     * bookkeeping and start clean.
     */
    public void resetSegments() {
        Arrays.fill(this.segment, Integer.MAX_VALUE);
        this.segData.clear();
    }

    /**
     * Remember the currently processed continent's segment number.
     *
     * @param globalX X origin of collision on world map
     * @param globalY Y origin of collision on world map
     */
    public void selectCollisionSegment(int globalX, int globalY) {
        int index = this.getIndex(globalX, globalY);
        this.activeContinent = this.segData.size();
        if (DEBUG) {
            if (index >= this.width * this.height) {
                System.out.println("Collision segment cannot be set outside plate!");
                System.exit(1);
            }
        }
        this.activeContinent = this.segment[index];
        if (DEBUG) {
            if (this.activeContinent >= this.segData.size()) {
                System.out.println("Collision happened at unsegmented location!");
                System.exit(1);
            }
        }
    }

    /**
     * Set the amount of plate's crustal material at some location.
     * <p>
     * If amount of crust to be set is negative, it'll be set to zero.
     *
     * @param globalX X offset on world map
     * @param globalY y offset on world map
     * @param crust   Amount of crust at the given location
     * @param age     Time of creation of new crust
     */
    public void setCrust(int globalX, int globalY, float crust, int age) {
        // Do not accept negative values.
        if (crust < 0) {
            crust = 0;
        }
        int index = this.getIndex(globalX, globalY);
        if (index < 0 || index >= this.width * this.height) {
            System.out.println("Extending plate");
            if (DEBUG) {
                if (crust <= 0) {
                    System.out.println("Extending plate for nothing!");
                    System.exit(1);
                }
            }
            int left = (int) this.x;
            int top = (int) this.y;
            int right = left + this.width - 1;
            int bottom = top + this.height - 1;
            globalX &= this.worldSize - 1;
            globalY &= this.worldSize - 1;
            // Calculate distance of new point from plate edges.
            int left0 = left - globalX;
            int right0 = (this.worldSize & -(globalX < left ? 1 : 0)) + globalX - right;
            int top0 = top - globalY;
            int bottom0 = (this.worldSize & -(globalY < top ? 1 : 0)) + globalY - bottom;
            if (left0 < 0) {
                left0 += Short.MAX_VALUE;
            }
            if (right0 < 0) {
                right0 += Short.MAX_VALUE;
            }
            if (top0 < 0) {
                top0 += Short.MAX_VALUE;
            }
            if (bottom0 < 0) {
                bottom0 += Short.MAX_VALUE;
            }
            // Set larger of horizontal/vertical distance to zero.
            // A valid distance is NEVER larger than world's side's length!
            int d_lft = left0 & -(left0 < right0 ? 1 : 0) & -(left0 < this.worldSize ? 1 : 0);
            int d_rgt = right0 & -(right0 <= left0 ? 1 : 0) & -(right0 < this.worldSize ? 1 : 0);
            int d_top = top0 & -(top0 < bottom0 ? 1 : 0) & -(top0 < this.worldSize ? 1 : 0);
            int d_btm = bottom0 & -(bottom0 <= top0 ? 1 : 0) & -(bottom0 < this.worldSize ? 1 : 0);
            // Scale all changes to multiple of 8.
            d_lft = (d_lft > 0 ? 1 : 0) + (d_lft >> 3) << 3;
            d_rgt = (d_rgt > 0 ? 1 : 0) + (d_rgt >> 3) << 3;
            d_top = (d_top > 0 ? 1 : 0) + (d_top >> 3) << 3;
            d_btm = (d_btm > 0 ? 1 : 0) + (d_btm >> 3) << 3;
            // Make sure plate doesn't grow bigger than the system it's in!
            if (this.width + d_lft + d_rgt > this.worldSize) {
                d_lft = 0;
                d_rgt = this.worldSize - this.width;
            }
            if (this.height + d_top + d_btm > this.worldSize) {
                d_top = 0;
                d_btm = this.worldSize - this.height;
            }
            if (DEBUG) {
                if (d_lft + d_top + d_rgt + d_btm == 0) {
                    System.out.printf("[%d, %d]x[%d, %d], [%d, %d]/[%d, %d]\n", (int) this.x, (int) this.y, (int) this.x + this.width, (int) this.y + this.height, globalX + this.worldSize * (globalX < this.worldSize ? 1 : 0), globalY + this.worldSize * (globalY < this.worldSize ? 1 : 0), globalX % this.worldSize, globalY % this.worldSize);
                    System.out.println("Index out of bounds, but nowhere to grow!");
                    System.exit(1);
                }
            }
            final int old_width = this.width;
            final int old_height = this.height;
            this.x -= d_lft;
            this.x += this.x >= 0 ? 0 : this.worldSize;
            this.width += d_lft + d_rgt;
            this.y -= d_top;
            this.y += this.y >= 0 ? 0 : this.worldSize;
            this.height += d_top + d_btm;
            float[] tmph = new float[this.width * this.height];
            int[] tmpa = new int[this.width * this.height];
            int[] tmps = new int[this.width * this.height];
            Arrays.fill(tmps, Integer.MAX_VALUE);
            // copy old plate into new.
            for (int j = 0; j < old_height; ++j) {
                final int dest_i = (d_top + j) * this.width + d_lft;
                final int src_i = j * old_width;
                System.arraycopy(this.heightmap, src_i, tmph, dest_i, old_width);
                System.arraycopy(this.age, src_i, tmpa, dest_i, old_width);
                System.arraycopy(this.segment, src_i, tmps, dest_i, old_width);
            }
            this.heightmap = tmph;
            this.age = tmpa;
            this.segment = tmps;
            // Shift all segment data to match new coordinates.
            for (int s = 0; s < this.segData.size(); ++s) {
                this.segData.get(s).x0 += d_lft;
                this.segData.get(s).x1 += d_lft;
                this.segData.get(s).y0 += d_top;
                this.segData.get(s).y1 += d_top;
            }
            index = this.getIndex(globalX, globalY);
            if (DEBUG) {
                if (index < 0 || index >= this.width * this.height) {
                    System.out.printf("Index out of bounds after resize!\n[%d, %d]x[%d, %d], [%d, %d]/[%d, %d]\n", (int) this.x, (int) this.y, (int) this.x + this.width, (int) this.y + this.height, globalX, globalY, globalX % this.worldSize, globalY % this.worldSize);
                    System.exit(1);
                }
            }
        }
        // Update crust's age.
        // If old crust exists, new age is mean of original and supplied ages.
        // If no new crust is added, original time remains intact.
        int oldCrust = -(this.heightmap[index] > 0 ? 1 : 0);
        int newCrust = -(crust > 0 ? 1 : 0);
        age = age & ~oldCrust | (int) ((this.heightmap[index] * this.age[index] + crust * age) / (this.heightmap[index] + crust)) & oldCrust;
        this.age[index] = age & newCrust | this.age[index] & ~newCrust;
        this.mass -= this.heightmap[index];
        // Set new crust height to desired location.
        this.heightmap[index] = crust;
        // Update mass counter.
        this.mass += crust;
    }

    /**
     * Container for details about a segmented crust area on this plate.
     */
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
