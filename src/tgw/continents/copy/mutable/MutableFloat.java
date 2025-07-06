package tgw.continents.copy.mutable;

public class MutableFloat {

    private float value;

    public MutableFloat() {
    }

    public MutableFloat(float value) {
        this.value = value;
    }

    public float get() {
        return this.value;
    }

    public void set(float value) {
        this.value = value;
    }
}
