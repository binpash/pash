for i in '#'
do
        cat << EOF
        x=\`printf '%s' \\$i\`; printf '%s\\n' "\$x"
EOF
done
