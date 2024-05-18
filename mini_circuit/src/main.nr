fn main(x: Field, y: pub Field) {
    assert(x * 4 == y * 3);
    ligma(x, y);
}

fn ligma(x: Field, y: Field) {
    assert(y - x == 2);
}

#[test]
fn test_main() {
    main(6, 8);

    // Uncomment to make test fail
    // main(1, 1);
}
