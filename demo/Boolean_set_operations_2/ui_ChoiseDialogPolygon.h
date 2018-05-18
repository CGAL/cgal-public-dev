/********************************************************************************
** Form generated from reading UI file 'ChoiseDialogPolygon.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CHOISEDIALOGPOLYGON_H
#define UI_CHOISEDIALOGPOLYGON_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QListView>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_ChoiseDialogPolygon
{
public:
    QVBoxLayout *verticalLayout;
    QListView *listView;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *ChoiseDialogPolygon)
    {
        if (ChoiseDialogPolygon->objectName().isEmpty())
            ChoiseDialogPolygon->setObjectName(QStringLiteral("ChoiseDialogPolygon"));
        ChoiseDialogPolygon->resize(400, 300);
        verticalLayout = new QVBoxLayout(ChoiseDialogPolygon);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        listView = new QListView(ChoiseDialogPolygon);
        listView->setObjectName(QStringLiteral("listView"));
        listView->setEditTriggers(QAbstractItemView::NoEditTriggers);
        listView->setAlternatingRowColors(false);
        listView->setUniformItemSizes(true);

        verticalLayout->addWidget(listView);

        buttonBox = new QDialogButtonBox(ChoiseDialogPolygon);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(ChoiseDialogPolygon);
        QObject::connect(buttonBox, SIGNAL(accepted()), ChoiseDialogPolygon, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), ChoiseDialogPolygon, SLOT(reject()));

        QMetaObject::connectSlotsByName(ChoiseDialogPolygon);
    } // setupUi

    void retranslateUi(QDialog *ChoiseDialogPolygon)
    {
        ChoiseDialogPolygon->setWindowTitle(QApplication::translate("ChoiseDialogPolygon", "Dialog", 0));
    } // retranslateUi

};

namespace Ui {
    class ChoiseDialogPolygon: public Ui_ChoiseDialogPolygon {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CHOISEDIALOGPOLYGON_H
