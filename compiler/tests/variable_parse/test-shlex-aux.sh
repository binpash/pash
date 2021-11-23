comment_fun()
{
    cat > /dev/null #Consume data from pipe so writers don't get SIGPIPE
}

bad_quote_fun()
{
    echo ${asf"asd}
}
