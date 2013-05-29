for f in *.ply; do
    tail -n +11 $f >> ../sheep_full_centered.ply
done
