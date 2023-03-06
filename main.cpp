#include <iostream>
#include <QApplication>
#include "gui.h"

using namespace ui_gui;

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    gui guiWidget;
    guiWidget.show();

    return app.exec();
}
