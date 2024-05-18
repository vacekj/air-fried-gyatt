fn main(x: Field, y: pub Field, z: pub Field) {
    assert(x * 4 == y * 3);
    ligma(x, y, z);
}

fn ligma(x: Field, y: Field, z: Field) {
    assert(y - x + z == 2);
    let a = x * y;
    assert(a * 3 * z * x == 0);
}

#[test]
fn test_main() {
    main(6, 8, 0);

    // Uncomment to make test fail
    // main(1, 1, 0);
}
