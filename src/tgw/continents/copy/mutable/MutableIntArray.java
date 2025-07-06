package tgw.continents.copy.mutable;

import org.jetbrains.annotations.Nullable;

public class MutableIntArray {

    private int @Nullable [] array;

    public MutableIntArray() {
    }

    public MutableIntArray(int @Nullable [] array) {
        this.array = array;
    }

    public int @Nullable [] get() {
        return this.array;
    }

    public void set(int @Nullable [] array) {
        this.array = array;
    }
}
