package tgw.continents.copy.mutable;

import org.jetbrains.annotations.Nullable;

public class MutableFloatArray {

    private float @Nullable [] array;

    public MutableFloatArray(float @Nullable [] array) {
        this.array = array;
    }

    public MutableFloatArray() {
    }

    public float @Nullable [] get() {
        return this.array;
    }

    public void set(float @Nullable [] array) {
        this.array = array;
    }
}
