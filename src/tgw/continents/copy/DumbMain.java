package tgw.continents.copy;

import tgw.continents.FastRandom;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public final class DumbMain {

    private static final float COLOR_STEP = 1.5f;
    private static final int NUM_PLATES = 10;
    private static final int[] COLORS = new int[NUM_PLATES];
    public static final FastRandom RANDOM = new FastRandom("Mauren".hashCode());

    private DumbMain() {}

    public static void main(String[] args) throws IOException {
        for (int i = 0, len = COLORS.length; i < len; i++) {
            COLORS[i] = 0xFF00_0000 | RANDOM.nextInt(0xFF_FFFF);
        }
        Lithosphere lithosphere = new Lithosphere(4_096, 0.65f, 60, 0.001f, 5_000, 0.10f, 2);
        lithosphere.createPlates(NUM_PLATES);
        int updates = 0;
        save(lithosphere, updates);
        while (lithosphere.getPlateCount() != 0) {
            save(lithosphere, ++updates);
            lithosphere.update();
        }
        save(lithosphere, ++updates);
        System.out.println("Finished");
    }

    private static int pack(float r, float g, float b) {
        if (r < 0 || r > 1) {
            throw new IllegalArgumentException("Invalid red: " + r);
        }
        if (g < 0 || g > 1) {
            throw new IllegalArgumentException("Invalid red: " + g);
        }
        if (b < 0 || b > 1) {
            throw new IllegalArgumentException("Invalid red: " + b);
        }
        return 0xFF << 24 | (int) (r * 255) << 16 | (int) (g * 255) << 8 | (int) (b * 255);
    }

    public static void save(Lithosphere lithosphere, int id) throws IOException {
        System.out.println("Saving " + id);
        float[] topography = lithosphere.getTopography();
        int sideLength = lithosphere.getSideLength();
        BufferedImage imageTopography = new BufferedImage(sideLength, sideLength, BufferedImage.TYPE_INT_ARGB);
        int i = 0;
        // Topographic color map:
        // 0.0:   0,   0,   0
        // 0.5:   0,   0, 255
        // 1.0:   0, 255, 255
        //
        // 0.0:   0, 255,   0
        // 1.0: 255, 255,   0
        // 2.5: 255, 128,   0
        // 3.0: 255,   0,   0
        // 4.0: 128,   0,   0
        // 5.0: 128,   0, 128 // Astral level 1
        // 6.0:   0,   0,   0 // Astral level 2
        for (int y = 0; y < sideLength; ++y) {
            for (int x = 0; x < sideLength; ++x) {
                float c = topography[i++];
                if (c < 0.5) {
                    imageTopography.setRGB(x, y, pack(0, 0, 2 * Math.max(0, c)));
                }
                else if (c < 1.0) {
                    imageTopography.setRGB(x, y, pack(0, 2 * (c - 0.5f), 1.0f));
                }
                else {
                    c -= 1.0f;
                    if (c < 1.0 * COLOR_STEP) {
                        imageTopography.setRGB(x, y, pack(c / COLOR_STEP, 1.0f, 0));
                    }
                    else if (c < 2.5 * COLOR_STEP) {
                        imageTopography.setRGB(x, y, pack(1.0f, 1.0f - 0.5f * (c - COLOR_STEP) / (1.5f * COLOR_STEP), 0));
                    }
                    else if (c < 3.0 * COLOR_STEP) {
                        imageTopography.setRGB(x, y, pack(1.0f, 0.5f - (c - 2.5f * COLOR_STEP) / COLOR_STEP, 0));
                    }
                    else if (c < 4.0 * COLOR_STEP) {
                        imageTopography.setRGB(x, y, pack(1.0f - 0.5f * (c - 3.0f * COLOR_STEP) / COLOR_STEP, 0, 0));
                    }
                    else if (c < 5.0 * COLOR_STEP) {
                        imageTopography.setRGB(x, y, pack(0.5f, 0, 0.5f * (c - 4.0f * COLOR_STEP) / COLOR_STEP));
                    }
                    else if (c < 6.0 * COLOR_STEP) {
                        imageTopography.setRGB(x, y, pack(0.5f - 0.5f * (c - 5.0f * COLOR_STEP) / COLOR_STEP, 0, 0.5f - 0.5f * (c - 5.0f * COLOR_STEP) / COLOR_STEP));
                    }
                    else {
                        imageTopography.setRGB(x, y, pack(0, 0, 0));
                    }
                }
            }
        }
        ImageIO.write(imageTopography, "png", new File(String.format("B:\\Desktop\\continents\\topography_%02d.png", id)));
        int[] plates = lithosphere.getPlates();
        BufferedImage imagePlates = new BufferedImage(sideLength, sideLength, BufferedImage.TYPE_INT_ARGB);
        i = 0;
        for (int y = 0; y < sideLength; ++y) {
            for (int x = 0; x < sideLength; ++x) {
                int c = plates[i++];
                imagePlates.setRGB(x, y, COLORS[c]);
            }
        }
        ImageIO.write(imagePlates, "png", new File(String.format("B:\\Desktop\\continents\\plates_%02d.png", id)));
    }
}
