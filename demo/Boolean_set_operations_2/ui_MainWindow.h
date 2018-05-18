/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QLocale>
#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QTreeView>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionSave;
    QAction *actionLoad;
    QAction *actionErase;
    QAction *actionDuplicate;
    QAction *actionEraseAll;
    QAction *actionSaveAs;
    QAction *actionShowHide;
    QAction *actionRecenterScene;
    QAction *actionSelectAllItems;
    QAction *actionDraw_Bezier_curve;
    QAction *actionDraw_Circular_curve;
    QAction *actionComplement;
    QAction *actionintersection;
    QAction *actionunion;
    QAction *actionDifference;
    QAction *actionSymmetric_Difference;
    QAction *actionMinkowski_Sum;
    QAction *actionHelp;
    QAction *actionDraw_Polygon;
    QAction *action_Quit;
    QWidget *centralwidget;
    QVBoxLayout *verticalLayout_2;
    QGraphicsView *SceneGraphicsView;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menuOperations;
    QMenu *menuView;
    QMenu *menu_Help;
    QStatusBar *statusbar;
    QDockWidget *sceneDockWidget;
    QWidget *dockWidgetContent;
    QGridLayout *gridLayout_2;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QToolButton *addButton;
    QToolButton *removeButton;
    QToolButton *duplicateButton;
    QToolButton *upButton;
    QToolButton *downButton;
    QSpacerItem *horizontalSpacer;
    QLineEdit *searchEdit;
    QTreeView *sceneView;
    QDockWidget *consoleDockWidget;
    QWidget *dockWidgetContents;
    QVBoxLayout *verticalLayout_3;
    QTextEdit *consoleTextEdit;
    QDockWidget *infoDockWidget;
    QWidget *dockWidgetContents_2;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout_4;
    QScrollArea *scrollArea;
    QWidget *scrollAreaWidgetContents;
    QGridLayout *gridLayout;
    QLabel *infoLabel;
    QLabel *displayLabel;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(978, 594);
        QIcon icon;
        icon.addFile(QStringLiteral(":/cgal/icons/resources/cgal_logo.xpm"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        actionSave = new QAction(MainWindow);
        actionSave->setObjectName(QStringLiteral("actionSave"));
        actionLoad = new QAction(MainWindow);
        actionLoad->setObjectName(QStringLiteral("actionLoad"));
        QIcon icon1;
        icon1.addFile(QStringLiteral(":/cgal/icons/plus"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoad->setIcon(icon1);
        actionErase = new QAction(MainWindow);
        actionErase->setObjectName(QStringLiteral("actionErase"));
        QIcon icon2;
        icon2.addFile(QStringLiteral(":/cgal/icons/minus"), QSize(), QIcon::Normal, QIcon::Off);
        actionErase->setIcon(icon2);
        actionDuplicate = new QAction(MainWindow);
        actionDuplicate->setObjectName(QStringLiteral("actionDuplicate"));
        QIcon icon3;
        icon3.addFile(QStringLiteral(":/cgal/icons/duplicate"), QSize(), QIcon::Normal, QIcon::Off);
        actionDuplicate->setIcon(icon3);
        actionEraseAll = new QAction(MainWindow);
        actionEraseAll->setObjectName(QStringLiteral("actionEraseAll"));
        actionSaveAs = new QAction(MainWindow);
        actionSaveAs->setObjectName(QStringLiteral("actionSaveAs"));
        actionShowHide = new QAction(MainWindow);
        actionShowHide->setObjectName(QStringLiteral("actionShowHide"));
        actionRecenterScene = new QAction(MainWindow);
        actionRecenterScene->setObjectName(QStringLiteral("actionRecenterScene"));
        actionSelectAllItems = new QAction(MainWindow);
        actionSelectAllItems->setObjectName(QStringLiteral("actionSelectAllItems"));
        actionDraw_Bezier_curve = new QAction(MainWindow);
        actionDraw_Bezier_curve->setObjectName(QStringLiteral("actionDraw_Bezier_curve"));
        QIcon icon4;
        icon4.addFile(QStringLiteral(":/cgal/icons/resources/insert_bezier.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionDraw_Bezier_curve->setIcon(icon4);
        actionDraw_Circular_curve = new QAction(MainWindow);
        actionDraw_Circular_curve->setObjectName(QStringLiteral("actionDraw_Circular_curve"));
        QIcon icon5;
        icon5.addFile(QStringLiteral(":/cgal/icons/resources/insert_polygon.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionDraw_Circular_curve->setIcon(icon5);
        actionComplement = new QAction(MainWindow);
        actionComplement->setObjectName(QStringLiteral("actionComplement"));
        QIcon icon6;
        icon6.addFile(QStringLiteral(":/cgal/icons/resources/comp_P.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionComplement->setIcon(icon6);
        actionintersection = new QAction(MainWindow);
        actionintersection->setObjectName(QStringLiteral("actionintersection"));
        QIcon icon7;
        icon7.addFile(QStringLiteral(":/cgal/icons/resources/intersection.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionintersection->setIcon(icon7);
        actionunion = new QAction(MainWindow);
        actionunion->setObjectName(QStringLiteral("actionunion"));
        QIcon icon8;
        icon8.addFile(QStringLiteral(":/cgal/icons/resources/union.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionunion->setIcon(icon8);
        actionDifference = new QAction(MainWindow);
        actionDifference->setObjectName(QStringLiteral("actionDifference"));
        QIcon icon9;
        icon9.addFile(QStringLiteral(":/cgal/icons/resources/diff_PQ.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionDifference->setIcon(icon9);
        actionSymmetric_Difference = new QAction(MainWindow);
        actionSymmetric_Difference->setObjectName(QStringLiteral("actionSymmetric_Difference"));
        QIcon icon10;
        icon10.addFile(QStringLiteral(":/cgal/icons/resources/symm_diff.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionSymmetric_Difference->setIcon(icon10);
        actionMinkowski_Sum = new QAction(MainWindow);
        actionMinkowski_Sum->setObjectName(QStringLiteral("actionMinkowski_Sum"));
        QIcon icon11;
        icon11.addFile(QStringLiteral(":/cgal/icons/resources/mink_sum.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionMinkowski_Sum->setIcon(icon11);
        actionHelp = new QAction(MainWindow);
        actionHelp->setObjectName(QStringLiteral("actionHelp"));
        actionDraw_Polygon = new QAction(MainWindow);
        actionDraw_Polygon->setObjectName(QStringLiteral("actionDraw_Polygon"));
        action_Quit = new QAction(MainWindow);
        action_Quit->setObjectName(QStringLiteral("action_Quit"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        verticalLayout_2 = new QVBoxLayout(centralwidget);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        SceneGraphicsView = new QGraphicsView(centralwidget);
        SceneGraphicsView->setObjectName(QStringLiteral("SceneGraphicsView"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(SceneGraphicsView->sizePolicy().hasHeightForWidth());
        SceneGraphicsView->setSizePolicy(sizePolicy);

        verticalLayout_2->addWidget(SceneGraphicsView);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 978, 22));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuEdit = new QMenu(menubar);
        menuEdit->setObjectName(QStringLiteral("menuEdit"));
        menuOperations = new QMenu(menubar);
        menuOperations->setObjectName(QStringLiteral("menuOperations"));
        menuView = new QMenu(menubar);
        menuView->setObjectName(QStringLiteral("menuView"));
        menu_Help = new QMenu(menubar);
        menu_Help->setObjectName(QStringLiteral("menu_Help"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        MainWindow->setStatusBar(statusbar);
        sceneDockWidget = new QDockWidget(MainWindow);
        sceneDockWidget->setObjectName(QStringLiteral("sceneDockWidget"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(sceneDockWidget->sizePolicy().hasHeightForWidth());
        sceneDockWidget->setSizePolicy(sizePolicy1);
        sceneDockWidget->setLocale(QLocale(QLocale::English, QLocale::UnitedStates));
        dockWidgetContent = new QWidget();
        dockWidgetContent->setObjectName(QStringLiteral("dockWidgetContent"));
        gridLayout_2 = new QGridLayout(dockWidgetContent);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        addButton = new QToolButton(dockWidgetContent);
        addButton->setObjectName(QStringLiteral("addButton"));
        addButton->setIcon(icon1);

        horizontalLayout->addWidget(addButton);

        removeButton = new QToolButton(dockWidgetContent);
        removeButton->setObjectName(QStringLiteral("removeButton"));
        removeButton->setIcon(icon2);

        horizontalLayout->addWidget(removeButton);

        duplicateButton = new QToolButton(dockWidgetContent);
        duplicateButton->setObjectName(QStringLiteral("duplicateButton"));
        duplicateButton->setIcon(icon3);

        horizontalLayout->addWidget(duplicateButton);

        upButton = new QToolButton(dockWidgetContent);
        upButton->setObjectName(QStringLiteral("upButton"));
        QIcon icon12;
        icon12.addFile(QStringLiteral(":/cgal/icons/resources/up.png"), QSize(), QIcon::Normal, QIcon::Off);
        upButton->setIcon(icon12);

        horizontalLayout->addWidget(upButton);

        downButton = new QToolButton(dockWidgetContent);
        downButton->setObjectName(QStringLiteral("downButton"));
        QIcon icon13;
        icon13.addFile(QStringLiteral(":/cgal/icons/resources/down.png"), QSize(), QIcon::Normal, QIcon::Off);
        downButton->setIcon(icon13);

        horizontalLayout->addWidget(downButton);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        searchEdit = new QLineEdit(dockWidgetContent);
        searchEdit->setObjectName(QStringLiteral("searchEdit"));

        horizontalLayout->addWidget(searchEdit);


        verticalLayout->addLayout(horizontalLayout);

        sceneView = new QTreeView(dockWidgetContent);
        sceneView->setObjectName(QStringLiteral("sceneView"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(sceneView->sizePolicy().hasHeightForWidth());
        sceneView->setSizePolicy(sizePolicy2);
        sceneView->setEditTriggers(QAbstractItemView::DoubleClicked|QAbstractItemView::EditKeyPressed|QAbstractItemView::SelectedClicked);
        sceneView->setDragDropMode(QAbstractItemView::InternalMove);
        sceneView->setAlternatingRowColors(true);
        sceneView->setSelectionMode(QAbstractItemView::ExtendedSelection);
        sceneView->setSelectionBehavior(QAbstractItemView::SelectRows);
        sceneView->setIndentation(10);

        verticalLayout->addWidget(sceneView);


        gridLayout_2->addLayout(verticalLayout, 0, 0, 1, 1);

        sceneDockWidget->setWidget(dockWidgetContent);
        MainWindow->addDockWidget(static_cast<Qt::DockWidgetArea>(1), sceneDockWidget);
        consoleDockWidget = new QDockWidget(MainWindow);
        consoleDockWidget->setObjectName(QStringLiteral("consoleDockWidget"));
        consoleDockWidget->setAllowedAreas(Qt::BottomDockWidgetArea|Qt::LeftDockWidgetArea|Qt::TopDockWidgetArea);
        dockWidgetContents = new QWidget();
        dockWidgetContents->setObjectName(QStringLiteral("dockWidgetContents"));
        verticalLayout_3 = new QVBoxLayout(dockWidgetContents);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        consoleTextEdit = new QTextEdit(dockWidgetContents);
        consoleTextEdit->setObjectName(QStringLiteral("consoleTextEdit"));
        consoleTextEdit->setReadOnly(true);

        verticalLayout_3->addWidget(consoleTextEdit);

        consoleDockWidget->setWidget(dockWidgetContents);
        MainWindow->addDockWidget(static_cast<Qt::DockWidgetArea>(8), consoleDockWidget);
        infoDockWidget = new QDockWidget(MainWindow);
        infoDockWidget->setObjectName(QStringLiteral("infoDockWidget"));
        dockWidgetContents_2 = new QWidget();
        dockWidgetContents_2->setObjectName(QStringLiteral("dockWidgetContents_2"));
        horizontalLayout_2 = new QHBoxLayout(dockWidgetContents_2);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(0);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        scrollArea = new QScrollArea(dockWidgetContents_2);
        scrollArea->setObjectName(QStringLiteral("scrollArea"));
        scrollArea->setMinimumSize(QSize(350, 0));
        scrollArea->setFrameShape(QFrame::NoFrame);
        scrollArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QStringLiteral("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 541, 173));
        gridLayout = new QGridLayout(scrollAreaWidgetContents);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        infoLabel = new QLabel(scrollAreaWidgetContents);
        infoLabel->setObjectName(QStringLiteral("infoLabel"));
        sizePolicy.setHeightForWidth(infoLabel->sizePolicy().hasHeightForWidth());
        infoLabel->setSizePolicy(sizePolicy);
        infoLabel->setContextMenuPolicy(Qt::DefaultContextMenu);
        infoLabel->setLineWidth(0);
        infoLabel->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        infoLabel->setTextInteractionFlags(Qt::LinksAccessibleByKeyboard|Qt::LinksAccessibleByMouse|Qt::TextBrowserInteraction|Qt::TextSelectableByKeyboard|Qt::TextSelectableByMouse);

        gridLayout->addWidget(infoLabel, 0, 0, 1, 1);

        scrollArea->setWidget(scrollAreaWidgetContents);

        verticalLayout_4->addWidget(scrollArea);

        displayLabel = new QLabel(dockWidgetContents_2);
        displayLabel->setObjectName(QStringLiteral("displayLabel"));

        verticalLayout_4->addWidget(displayLabel);


        horizontalLayout_2->addLayout(verticalLayout_4);

        infoDockWidget->setWidget(dockWidgetContents_2);
        MainWindow->addDockWidget(static_cast<Qt::DockWidgetArea>(8), infoDockWidget);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuEdit->menuAction());
        menubar->addAction(menuOperations->menuAction());
        menubar->addAction(menuView->menuAction());
        menubar->addAction(menu_Help->menuAction());
        menuFile->addAction(actionLoad);
        menuFile->addAction(actionErase);
        menuFile->addAction(actionEraseAll);
        menuFile->addAction(actionDuplicate);
        menuFile->addAction(actionSaveAs);
        menuFile->addAction(actionSave);
        menuFile->addSeparator();
        menuFile->addAction(action_Quit);
        menuEdit->addAction(actionShowHide);
        menuEdit->addAction(actionSelectAllItems);
        menuOperations->addAction(actionDraw_Bezier_curve);
        menuOperations->addAction(actionDraw_Circular_curve);
        menuOperations->addAction(actionComplement);
        menuOperations->addAction(actionintersection);
        menuOperations->addAction(actionunion);
        menuOperations->addAction(actionDifference);
        menuOperations->addAction(actionSymmetric_Difference);
        menuOperations->addAction(actionMinkowski_Sum);
        menuView->addAction(actionRecenterScene);
        menu_Help->addAction(actionHelp);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Boolean Set Operations and Minkowksi Sum", 0));
        actionSave->setText(QApplication::translate("MainWindow", "Save", 0));
        actionSave->setShortcut(QApplication::translate("MainWindow", "Ctrl+Q", 0));
        actionLoad->setText(QApplication::translate("MainWindow", "&Load curve", 0));
        actionLoad->setShortcut(QApplication::translate("MainWindow", "Ctrl+L", 0));
        actionErase->setText(QApplication::translate("MainWindow", "&Erase curve", 0));
        actionErase->setShortcut(QApplication::translate("MainWindow", "Del", 0));
        actionDuplicate->setText(QApplication::translate("MainWindow", "&Duplicate curve", 0));
        actionDuplicate->setShortcut(QApplication::translate("MainWindow", "Ctrl+D", 0));
        actionEraseAll->setText(QApplication::translate("MainWindow", "&Erase All curves", 0));
        actionSaveAs->setText(QApplication::translate("MainWindow", "Save &as...", 0));
        actionShowHide->setText(QApplication::translate("MainWindow", "Show/Hide Curves", 0));
        actionShowHide->setShortcut(QApplication::translate("MainWindow", "Ctrl+Space", 0));
        actionRecenterScene->setText(QApplication::translate("MainWindow", "Re&center Scene", 0));
        actionRecenterScene->setShortcut(QApplication::translate("MainWindow", "Ctrl+R", 0));
        actionSelectAllItems->setText(QApplication::translate("MainWindow", "Select All Curves", 0));
        actionSelectAllItems->setShortcut(QApplication::translate("MainWindow", "Ctrl+A", 0));
        actionDraw_Bezier_curve->setText(QApplication::translate("MainWindow", "Draw Bezier polygon", 0));
        actionDraw_Circular_curve->setText(QApplication::translate("MainWindow", "Draw Circular polygon", 0));
        actionComplement->setText(QApplication::translate("MainWindow", "Complement", 0));
        actionintersection->setText(QApplication::translate("MainWindow", "Intersection", 0));
        actionunion->setText(QApplication::translate("MainWindow", "Union", 0));
        actionDifference->setText(QApplication::translate("MainWindow", "Difference", 0));
        actionSymmetric_Difference->setText(QApplication::translate("MainWindow", "Symmetric Difference", 0));
        actionMinkowski_Sum->setText(QApplication::translate("MainWindow", "Minkowski Sum", 0));
        actionHelp->setText(QApplication::translate("MainWindow", "&About", 0));
        actionDraw_Polygon->setText(QApplication::translate("MainWindow", "Draw Polygon", 0));
        action_Quit->setText(QApplication::translate("MainWindow", "&Quit", 0));
        menuFile->setTitle(QApplication::translate("MainWindow", "&File", 0));
        menuEdit->setTitle(QApplication::translate("MainWindow", "&Edit", 0));
        menuOperations->setTitle(QApplication::translate("MainWindow", "&Operations", 0));
        menuView->setTitle(QApplication::translate("MainWindow", "&View", 0));
        menu_Help->setTitle(QApplication::translate("MainWindow", "&Help", 0));
        sceneDockWidget->setWindowTitle(QApplication::translate("MainWindow", "Geometric Objects", 0));
        addButton->setText(QApplication::translate("MainWindow", "+", 0));
        removeButton->setText(QApplication::translate("MainWindow", "-", 0));
        duplicateButton->setText(QApplication::translate("MainWindow", "...", 0));
        upButton->setText(QApplication::translate("MainWindow", "...", 0));
        downButton->setText(QApplication::translate("MainWindow", "...", 0));
        searchEdit->setText(QApplication::translate("MainWindow", "no Idea", 0));
        searchEdit->setPlaceholderText(QApplication::translate("MainWindow", "Quick filter... <Alt+Q>", 0));
        consoleDockWidget->setWindowTitle(QApplication::translate("MainWindow", "Console", 0));
        infoDockWidget->setWindowTitle(QApplication::translate("MainWindow", "Infos", 0));
        displayLabel->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
