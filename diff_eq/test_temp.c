#include <stdlib.h>
#include <unistd.h>

int main()
{
    system("cat /sys/class/thermal/thermal_zone0/temp > tempGraph");
    for(;;)
    {
        system("cat /sys/class/thermal/thermal_zone0/temp >> tempGraph");
        sleep(5);
    }
}