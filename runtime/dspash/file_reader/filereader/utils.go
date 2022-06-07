package filereader

import (
	"fmt"
	"os/exec"
)

func GetAbsPath(s string) (string, error) {
	// Hacky but should work for now
	out, err := exec.Command("bash", "-c", fmt.Sprintf("echo -n %s", s)).Output()
	if err != nil {
		return "", err
	}
	return string(out), nil
}
