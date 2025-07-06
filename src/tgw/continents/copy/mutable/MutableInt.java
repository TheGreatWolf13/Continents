package tgw.continents.copy.mutable;

public class MutableInt {

    private int value;

    public MutableInt(int value) {
        this.value = value;
    }

    public MutableInt() {
    }

    public int get() {
        return this.value;
    }

    public void set(int value) {
        this.value = value;
    }
}
