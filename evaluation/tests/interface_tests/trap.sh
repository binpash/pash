myfunction()
{
    echo myfunction invoked
}
trap myfunction EXIT
echo hello one
echo hello two
